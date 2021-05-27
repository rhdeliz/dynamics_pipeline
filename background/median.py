import time
import numpy as np
from scipy import ndimage
import tifffile, ray
import settings

@ray.remote
def median_blur(img, median_filter_size):
    med = ndimage.median_filter(img, size=median_filter_size)
    med_rm = img - np.minimum(med, img)

    return med_rm

def remove_median_blur(protein_image_x, background_remove_image_list, cell_size, puncta_size, tiff_compression_level):
    try:
        # Get protein info
        original_image_path = background_remove_image_list['df_rm_save_path'][protein_image_x]
        cell_med_rm_save_path = background_remove_image_list['cell_med_rm_save_path'][protein_image_x]

        # Import images
        original_image = tifffile.imread(original_image_path)

        # Get frame count
        frames = original_image.shape[0]
        frames = range(0, frames)

        # Run in parallel
        ray.init()
        result_ids = []
        for frame in frames:
            result_id = median_blur.remote(original_image[frame], cell_size)
            result_ids.append(result_id)
        med_rm_img = settings.parallel.ids_to_vals(result_ids)
        # Save image
        tifffile.imsave(cell_med_rm_save_path, med_rm_img, bigtiff=True, compress=tiff_compression_level, dtype=med_rm_img[0].dtype)
        time.sleep(5)
        ray.shutdown()

        # Import images
        original_image = tifffile.imread(cell_med_rm_save_path)
        puncta_med_rm_save_path = background_remove_image_list['puncta_med_rm_save_path'][protein_image_x]

        # Run in parallel
        ray.init()
        result_ids = []
        for frame in frames:
            result_id = median_blur.remote(original_image[frame], puncta_size)
            result_ids.append(result_id)
        med_rm_img = settings.parallel.ids_to_vals(result_ids)
        # Save image
        tifffile.imsave(puncta_med_rm_save_path, med_rm_img, bigtiff=True, compress=tiff_compression_level, dtype=med_rm_img[0].dtype)
        time.sleep(5)
        ray.shutdown()

    except:
        print("Cannot complete remove")
    return original_image_path
