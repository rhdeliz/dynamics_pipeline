'''
Converts ND2 file to multi-page TIFF
'''

import shutil, warnings
from pims_nd2 import ND2_Reader  # Install as pims-nd2
import ray, tifffile, os, time
import numpy as np
import cv2 as cv
from scipy import ndimage
from subprocess import call


@ray.remote
def make_tiff(image_x_path, img_processing_path, tiff_path, ligand_puncta_radius):
    try:
        warnings.filterwarnings("ignore")
        img = ND2_Reader(image_x_path)
        warnings.filterwarnings("default")
    except:
        print('Cannot open ND2 image')

    try:
        save_path = img_processing_path
        # Create if it doesn't exist
        if not os.path.exists(save_path):
            os.makedirs(save_path)
    except:
        print('Cannot create image path')

    try:
        # Get frame count
        n_frames = img.sizes['t'] - 1
        frames = range(0, n_frames)

        img_frames = []
        for t in frames:
            orig = img.get_frame_2D(c=0, x=0, y=0, t=t)
            img_frames.append(orig)

        # Avg pixel intensity
        img_min = np.min(img_frames, axis=0)
        img_subtracted = []
        for t in frames:
            img_frame = img_frames[t]
            img_frame = cv.subtract(img_frame, img_min)
            img_subtracted.append(img_frame)

        # Max pixel intensity
        img_max = np.max(img_subtracted, axis=0)
        kernel = cv.getStructuringElement(cv.MORPH_ELLIPSE, (ligand_puncta_radius*3, ligand_puncta_radius*3))
        dilated_img = cv.dilate(img_max, kernel, iterations=1)

        # Median blur
        median_img = ndimage.median_filter(dilated_img, size=ligand_puncta_radius*10)

        # Flat-field correct
        ff_img = []
        for t in frames:
            img_frame = img_subtracted[t]
            img_frame = img_frame*500
            img_frame = np.divide(img_frame, median_img)
            ff_img.append(img_frame)

        # Prepare for 8-bit
        max_value = np.max(ff_img)
        # Flat-field correct
        rescaled_img = []
        for t in frames:
            ff_img_frame = ff_img[t]
            ff_img_frame = ff_img_frame/max_value
            ff_img_frame = ff_img_frame*2**8
            rescaled_img.append(ff_img_frame)

        # Save file
        tifffile.imsave(tiff_path, rescaled_img, bigtiff=False, dtype='uint8')

    except:
        print('Cannot process image')

    try:
        # Move ND2 file once finished
        shutil.move(image_x_path, img_processing_path)
        # Clos e image
        img.close()
        time.sleep(5)
    except:
        print('Cannot move processed nd2 file')

    return image_x_path


def trackmate(imagej, protein_path, image_path, trackmate_threshold, trackmate_frame_gap, trackmate_max_link_distance, trackmate_gap_link_distance, puncta_diameter):
    # Get macro path
    macro_path = os.getcwd()
    macro_path = os.path.join(macro_path, 'ligand', 'run_trackmate.ijm')
    # Numbers to string
    trackmate_threshold = str(trackmate_threshold)
    trackmate_frame_gap = str(trackmate_frame_gap)
    trackmate_max_link_distance = str(trackmate_max_link_distance)
    trackmate_gap_link_distance = str(trackmate_gap_link_distance)
    puncta_diameter = str(puncta_diameter)

    # Execute
    separator = '\n'
    call([imagej, '--headless', '--console', '-macro', macro_path,
          protein_path + separator + image_path + separator +
          trackmate_threshold + separator + trackmate_frame_gap + separator + trackmate_max_link_distance + separator + trackmate_gap_link_distance + separator +
          puncta_diameter])


    return image_path
