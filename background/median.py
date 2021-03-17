import os, tifffile, math, time
import numpy as np
from pims_nd2 import ND2_Reader
from PIL import Image
from scipy import ndimage
from matplotlib import pyplot as plt

import cv2 as cv
from matplotlib import pyplot as plt

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
