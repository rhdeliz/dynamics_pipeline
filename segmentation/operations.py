import math, ray, tifffile
import numpy as np
import cv2 as cv
import os
import matplotlib.pyplot as plt
import time

@ray.remote
def segment(segmentation_input_path, cell_size, tiff_compression_level):
    # Import images
    img = tifffile.imread(segmentation_input_path)

    # Max pixel intensity
    img_max = np.max(img, axis=0)
    # Convert to 8-bit
    img_max = img_max - img_max.min()
    img_max = img_max / img_max.max() * 255
    img_max = img_max.astype(np.uint8)

    # Average pixel intensity
    img = np.mean(img, axis=0)
    img = np.round(img)
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
    blur_img = cv.blur(contrast_img, (cell_size*3, cell_size*3))
    blur_img = blur_img * 0.5
    blur_img = blur_img.astype(np.uint8)
    # Median image
    median_img = cv.medianBlur(contrast_img, ksize=cell_size)
    # Subtract blurs
    subtracted_img = cv.subtract(median_img, blur_img)

    # Dilate maxima
    kernel = cv.getStructuringElement(cv.MORPH_ELLIPSE, (cell_size*3, cell_size*3))
    dilated_img = cv.dilate(subtracted_img, kernel, iterations=1)

    # Gamma correct image again to make dim objects brighter
    mid = 1
    mean = np.mean(dilated_img)
    gamma = math.log(mid * 255) / math.log(mean)
    gamma_img = np.power(dilated_img, gamma).clip(0, 255).astype(np.uint8)

    # Make binary applying a threshold
    ret, mask_img = cv.threshold(gamma_img, 1, 255, cv.THRESH_BINARY + cv.THRESH_OTSU)
    mask_img_color = cv.cvtColor(mask_img, cv.COLOR_GRAY2BGR)

    # Create marker
    # Adjust contrast
    contrast_img = img_max.astype(np.uint16)
    contrast_img = contrast_img * 6.67
    contrast_img = np.clip(contrast_img, 0, 255)
    contrast_img = contrast_img - contrast_img.min()
    contrast_img = contrast_img / contrast_img.max() * 255
    contrast_img = contrast_img.astype(np.uint8)

    # Blur image
    blur_img = cv.blur(contrast_img, (cell_size*3, cell_size*3))
    # Median image
    median_img = cv.medianBlur(contrast_img, ksize=cell_size)
    # Subtract blurs
    subtracted_img = cv.subtract(median_img, blur_img)
    # Gamma correct image
    gamma_img = np.power(subtracted_img, 1.5).clip(0,255).astype(np.uint8)
    # Median blur to remove noise
    median_img_2 = cv.medianBlur(gamma_img, ksize=cell_size)

    # Marker labelling
    ret, marker_img = cv.connectedComponents(median_img_2)

    # Segment
    final_segmentation = cv.watershed(mask_img_color, marker_img)
    final_segmentation = final_segmentation.astype(np.uint8)
    # Eliminate background
    final_segmentation = marker_img * (mask_img / 255)
    final_segmentation = final_segmentation.clip(0,254)
    final_segmentation = final_segmentation.astype('uint8')
    # Save
    # Create folders
    image_path = os.path.dirname(segmentation_input_path)
    segmentation_output_path = os.path.join(image_path, 'Segmentation.tif')
    tifffile.imsave(segmentation_output_path, final_segmentation, bigtiff=True, compress=tiff_compression_level, dtype=  final_segmentation[0].dtype)

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
    contrast_img = contrast_img - np.mean(contrast_img)
    contrast_img = np.clip(contrast_img, 0, 255)
    contrast_img = contrast_img - contrast_img.min()
    contrast_img = contrast_img / contrast_img.max() * 255
    final_segmentation_img = contrast_img.astype(np.uint8)
    # Add boundaries
    final_segmentation_img[eroded_img == 255] = [255]

    plt.imshow(final_segmentation_img, cmap='inferno')
    overlay_file = os.path.join(image_path, 'Segmentation.png')
    plt.savefig(overlay_file)

def process_cell_frame(frame_x, input_image, cell_x_mask, indexes):
    cell_x_img = input_image[frame_x]
    cutout_img = cell_x_img * cell_x_mask
    cutout_img = cutout_img[tuple(indexes)]

    return cutout_img

def process_cell(cell_x, segmentation_image, frames, input_image, image_x_path, area_results, tiff_compression_level):
    # Get cell masks
    cell_x_mask = segmentation_image == cell_x

    i, j = np.where(cell_x_mask)
    indexes = np.meshgrid(np.arange(min(i), max(i) + 1), np.arange(min(j), max(j) + 1), indexing='ij')

    # Run in parallel
    cropped_img = []
    for frame_x in frames:
        cropped_frame_img = process_cell_frame(frame_x, input_image, cell_x_mask, indexes)
        cropped_img.append(cropped_frame_img)
    # Save image
    save_path = os.path.dirname(image_x_path)
    cell_path = np.where(area_results == cell_x)[0][0]+1
    cell_path = 'Cell_' + cell_path.__str__()
    cell_path = os.path.join(save_path, cell_path)
    # Create cell path if it doesn't exist
    if not os.path.exists(cell_path):
        os.mkdir(cell_path)
    # Save name
    save_name = os.path.basename(image_x_path)
    save_path = os.path.join(cell_path, save_name)
    # Save image
    tifffile.imsave(save_path, cropped_img, bigtiff=True, compress=tiff_compression_level, dtype=cropped_img[0].dtype)
    time.sleep(5)
    return(save_name)

@ray.remote
def make_substacks(image_x_path, segmentation_image_x_path, puncta_size, tiff_compression_level):
    # Import images
    input_image = tifffile.imread(image_x_path)
    if os.path.exists(segmentation_image_x_path):
        segmentation_image = tifffile.imread(segmentation_image_x_path)
    else:
        segmentation_image = np.zeros((input_image.shape[1], input_image.shape[2]))+1

    # Get cell count
    n_cells = segmentation_image.max()
    n_cells = range(1, n_cells)

    # Get frame count
    frames = input_image.shape[0]
    frames = range(0, frames)

    # Get cells larger than puncta area
    area_results = []
    for cell_x in n_cells:
        result = sum(sum(segmentation_image == cell_x))
        result = result > puncta_size ** 2
        area_results.append(result)
    area_results = np.where(area_results)[0] + 1

    for cell_x in area_results:
        process_cell(cell_x, segmentation_image, frames, input_image, image_x_path, area_results, tiff_compression_level)

    time.sleep(5)
    ray.shutdown()
