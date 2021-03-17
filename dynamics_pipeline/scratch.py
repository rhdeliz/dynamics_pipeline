import os, tifffile, math, time
import numpy as np
from pims_nd2 import ND2_Reader
from PIL import Image
from scipy import ndimage
from matplotlib import pyplot as plt

import cv2 as cv
from matplotlib import pyplot as plt

def process_image(file):
    try:
        # Darkframe list
        df_path = '/Users/u_deliz/Desktop/NewPipeline/DarkFrames'
        df_Cy5 = 'AVG_20201026 Darkfield 50ms binned.tif'
        df_GFP = 'AVG_20201026 Darkfield 100ms binned 002.tif'
        df_RFP = 'AVG_20201026 Darkfield 200ms binned 003.tif'

        # Median filter size
        median_filter_size = 25

        # Open darkframe images
        try:
            os.chdir(df_path)
        except:
            print('df_path does not exist')

        try:
            df_Cy5 = Image.open(df_Cy5)
            df_Cy5 = np.array(df_Cy5)
        except:
            print('df_Cy5 does not exist')

        try:
            df_GFP = Image.open(df_GFP)
            df_GFP = np.array(df_GFP)
        except:
            print('df_GFP does not exist')

        try:
            df_RFP = Image.open(df_RFP)
            df_RFP = np.array(df_RFP)
        except:
            print('df_RFP does not exist')

        try:
            # Get image name and path
            img_path = os.path.abspath(os.path.join(file, os.pardir))

            # Image name
            save_name = os.path.basename(file)
            save_name_len = len(save_name) - 4
            save_name = save_name[:save_name_len]

            # Open image
            img = ND2_Reader(file)

            # Create new folder with image name
            img_dir = img_path + '/' + save_name
            try:
                os.mkdir(img_dir)
            except:
                pass
            os.chdir(img_dir)

            # Save metadata to text file
            try:
                file = open('Metadata.txt', 'w')
                file.write(img.metadata_text)
                file.close()
            except:
                print('Metadata error')

            # Function to get channel names
            def channel_names(img, planes):
                try:
                    channels = []
                    i = 0
                    while i <= planes:
                        plane_x = 'plane_' + str(i)
                        name = img.metadata[plane_x]['name']
                        channels.append(name)
                        i += 1
                    return channels
                except:
                    print('channel_names error')

            # Get channel names
            nChannels = img.metadata['plane_count'] - 1
            channels = channel_names(nChannels)
            nFrames = img.sizes['t'] - 1

            for channel in channels:
                #Function to get original
                def original(img, frames, c):
                    try:
                        t = 0
                        all_orig = []
                        while t <= frames:
                            orig = img.get_frame_2D(c=c, x=0, y=0, t=t)
                            orig = np.array(orig)
                            orig = Image.fromarray(orig)
                            all_orig.append(np.array(orig))
                            t += 1
                        return np.array(all_orig)
                    except:
                        print('original error')

                # Function to remove darkframe
                def remove_darkframe(img, frames, c, df_img):
                    try:
                        t = 0
                        all_df_rm = []
                        while t <= frames:
                            df_rm = img.get_frame_2D(c=c, x=0, y=0, t=t)
                            df_rm = np.array(df_rm)
                            df_rm = df_rm - np.minimum(df_rm, df_img)
                            #df_rm = cv.subtract(df_rm, df_img)
                            df_rm = Image.fromarray(df_rm)
                            all_df_rm.append(np.array(df_rm))
                            t += 1
                        return np.array(all_df_rm)
                    except:
                        print('remove_darkframe error')

                # Function to remove median
                def remove_median(img, frames):
                    try:
                        t = 0
                        all_med_rm = []
                        while t <= frames:
                            df_rm = img[t]
                            med = median_blur(df_rm, median_filter_size)
                            #med = ndimage.median_filter(df_rm, size=median_filter_size)
                            #med_rm = df_rm - np.minimum(med, df_rm)
                            med_rm = cv.subtract(df_rm, med)
                            med_rm = Image.fromarray(med_rm)
                            all_med_rm.append(np.array(med_rm))
                            t += 1
                        return np.array(all_med_rm)
                    except:
                        print('remove_median error')

                if channel == 'T Cy5':
                    try:
                        # Original
                        all_orig = original(img, nFrames, channels.index(channel))
                        out_name = save_name + ' C=' + str(channels.index(channel)) + '.tif'
                        tifffile.imsave(out_name, all_orig, compress=9)

                        # Remove darkframe
                        all_df_rm = remove_darkframe(img, nFrames, channels.index(channel), df_Cy5)
                        out_name = save_name + ' C=' + str(channels.index(channel)) + '-DFRm.tif'
                        tifffile.imsave(out_name, all_df_rm, compress=9)
                    except:
                        print('Error with ' + save_name + ' ' + channel)

                if channel == 'T GFP':
                    try:
                        # Original
                        all_orig = original(img, nFrames, channels.index(channel))
                        out_name = save_name + ' C=' + str(channels.index(channel)) + '.tif'
                        tifffile.imsave(out_name, all_orig, compress=9)

                        # Remove darkframe
                        all_df_rm = remove_darkframe(img, nFrames, channels.index(channel), df_GFP)
                        out_name = save_name + ' C=' + str(channels.index(channel)) + '-DFRm.tif'
                        tifffile.imsave(out_name, all_df_rm, compress=9)

                        # Remove median
                        all_med_rm = remove_median(all_df_rm, nFrames)
                        out_name = save_name + ' C=' + str(channels.index(channel)) + '-MedRm.tif'
                        tifffile.imsave(out_name, all_med_rm, compress=9)
                    except:
                        print('Error with ' + save_name + ' ' + channel)
                if channel == 'T RFP Cy3':
                    try:
                        # Original
                        all_orig = original(img, nFrames, channels.index(channel))
                        out_name = save_name + ' C=' + str(channels.index(channel)) + '.tif'
                        tifffile.imsave(out_name, all_orig, compress=9)

                        # Remove darkframe
                        all_df_rm = remove_darkframe(img, nFrames, channels.index(channel), df_RFP)
                        out_name = save_name + ' C=' + str(channels.index(channel)) + '-DFRm.tif'
                        tifffile.imsave(out_name, all_df_rm, compress=9)

                        # Remove median
                        all_med_rm = remove_median(all_df_rm, nFrames)
                        out_name = save_name + ' C=' + str(channels.index(channel)) + '-MedRm.tif'
                        tifffile.imsave(out_name, all_med_rm, compress=9)
                    except:
                        print('Error with ' + save_name + ' ' + channel)
        except:
            print('process_image error with ' + file)

    except:
        print('process_image error' + os.path.basename(file))


tic = time.time()
process_image('/Users/u_deliz/Desktop/NewPipeline/Images/20190705 GFP calibration.nd2')
toc = time.time()
print(toc - tic)


# FramesSequenceND.__init__()

def process_segmentation(img_file):
    # Import image
    img = cv.imreadmulti(img_file, flags=-1)
    # Max pixel intensity
    img = np.max(img[1], axis=0)
    # Convert to 8-bit
    img = img - img.min()
    img = img / img.max() * 255
    img = img.astype(np.uint8)

    # Create mask (background vs foreground)
    # Adjust contrast
    contrast_img = img.astype(np.uint16)
    contrast_img = contrast_img * 10
    contrast_img = np.clip(contrast_img, 0, 255)
    contrast_img = contrast_img - contrast_img.min()
    contrast_img = contrast_img / contrast_img.max() * 255
    contrast_img = contrast_img.astype(np.uint8)
    # Blur image
    blur_img = cv.blur(contrast_img, (75, 75))
    blur_img = blur_img * 0.5
    blur_img = blur_img.astype(np.uint8)
    # Median image
    median_img = cv.medianBlur(contrast_img, ksize=21)
    # Subtract blurs
    subtracted_img = cv.subtract(median_img, blur_img)
    # Dilate maxima
    kernel = cv.getStructuringElement(cv.MORPH_ELLIPSE, (75, 75))
    dilated_img = cv.dilate(subtracted_img, kernel, iterations=1)
    # Gamma correct image again to make dim objects brighter
    mid = 1.25
    mean = np.mean(dilated_img)
    gamma = math.log(mid * 255) / math.log(mean)
    gamma_img = np.power(dilated_img, gamma).clip(0, 255).astype(np.uint8)
    # Make binary applying a threshold
    ret, mask_img = cv.threshold(gamma_img, 1, 255, cv.THRESH_BINARY + cv.THRESH_OTSU)
    mask_img_color = cv.cvtColor(mask_img, cv.COLOR_GRAY2BGR)

    # Create marker
    # Adjust contrast
    contrast_img = img.astype(np.uint16)
    contrast_img = contrast_img * 5
    contrast_img = np.clip(contrast_img, 0, 255)
    contrast_img = contrast_img - contrast_img.min()
    contrast_img = contrast_img / contrast_img.max() * 255
    contrast_img = contrast_img.astype(np.uint8)
    # Blur image
    blur_img = cv.blur(contrast_img, (75, 75))
    # Median image
    median_img = cv.medianBlur(contrast_img, ksize=21)
    # Subtract blurs
    subtracted_img = cv.subtract(median_img, blur_img)
    # Gamma correct image
    mid = 0.05
    mean = np.mean(subtracted_img)
    gamma = math.log(mid * 255) / math.log(mean)
    gamma_img = np.power(subtracted_img, gamma).clip(0, 255).astype(np.uint8)
    # Adaptive threshold
    threshold_img = cv.adaptiveThreshold(gamma_img, 255, cv.ADAPTIVE_THRESH_MEAN_C, cv.THRESH_BINARY, 61, -10)
    # Median blur to remove noise
    median_img_2 = cv.medianBlur(threshold_img, ksize=21)
    # Marker labelling
    ret, marker_img = cv.connectedComponents(median_img_2)
    # Segment
    final_segmentation = cv.watershed(mask_img_color, marker_img)
    final_segmentation = final_segmentation.astype(np.uint8)
    # Eliminate background
    final_segmentation = marker_img * (mask_img / 255)

    # Plot over image
    # Find edges
    final_segmentation[final_segmentation > 0] = [255]
    # Erode Image
    kernel = cv.getStructuringElement(cv.MORPH_ELLIPSE, (3, 3))
    eroded_img = cv.erode(final_segmentation, kernel, iterations=1)
    eroded_img[eroded_img < 1] = [0]
    eroded_img = eroded_img.astype(np.uint8)
    # Subtract
    eroded_img = final_segmentation - eroded_img
    eroded_img[eroded_img != 0] = [255]
    # Adjust contrast
    contrast_img = img.astype(np.uint16)
    contrast_img = contrast_img * 3
    contrast_img = np.clip(contrast_img, 0, 255)
    contrast_img = contrast_img - contrast_img.min()
    contrast_img = contrast_img / contrast_img.max() * 255
    final_segmentation_img = contrast_img.astype(np.uint8)
    # Add boundaries
    final_segmentation_img[eroded_img == 255] = [255]
    # Plot
    plt.imshow(final_segmentation_img, cmap='inferno')

    # Save images
    save_folder
    Image.fromarray(final_segmentation_img)





