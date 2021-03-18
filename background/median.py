import time
import numpy as np
from scipy import ndimage
# from matplotlib import pyplot as plt
# from pims_nd2 import ND2_Reader
# from PIL import Image
# import os, tifffile, math,

import cv2 as cv


def remove_approxiamte_median_blur(img, median_filter_size):
    img2 = img / img.max() * 255
    img2 = img2.astype(np.uint8)
    med = cv.medianBlur(img2, ksize=median_filter_size)
    med = med.astype(np.uint16)
    med = med / 255 * img.max()
    med = np.round(med)
    med = med.astype(np.uint16)
    return med

def remove_real_median_blur(img, median_filter_size):
    med = ndimage.median_filter(img, size=median_filter_size)
    med_rm = img - np.minimum(med, img)

    return med_rm


import tifffile, ray
import cv2 as cv

@ray.remote
def remove(protein_image_x, background_remove_image_list, real_median, median_filter_size, tiff_compression_level):
    try:
        # Get protein info
        original_image_path = background_remove_image_list['df_rm_save_path'][protein_image_x]
        save_image_path = background_remove_image_list['med_rm_save_path'][protein_image_x]
        # Import images
        original_image = tifffile.imread(original_image_path)

        # Blur image
        # Get frame count
        frames = original_image.shape[0]
        frames = range(0, frames)

        if real_median == True:
            med = ndimage.median_filter(original_image, size=median_filter_size)

        else:
            rescaled_original_image = original_image / original_image.max() * 255
            rescaled_original_image = rescaled_original_image.astype(np.uint8)
            med = cv.medianBlur(rescaled_original_image, ksize=median_filter_size)
            med = med.astype(np.uint16)
            med = med / 255 * original_image.max()
            med = np.round(med)
            med = med.astype(np.uint16)




        # Subtract
        dark_frame_removed_image = []
        for t in frames:
            df_rm = original_image[t]
            # df_rm = df_rm - np.minimum(df_rm, dark_frame_image)
            df_rm = cv.subtract(df_rm, dark_frame_image)
            dark_frame_removed_image.append(df_rm)




        # Save image
        tifffile.imsave(save_image_path, dark_frame_removed_image, compress=tiff_compression_level)
        time.sleep(5)
    except:
        print("Cannot complete remove")
    return original_image_path

