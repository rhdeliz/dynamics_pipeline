import cv2 as cv # Install as opencv-python
import ray
import tifffile
import time
import background.parameters as parameters


# opencv-python
@ray.remote
def remove_frame(protein_image_x, background_remove_image_list, dark_frames_list, tiff_compression_level):
    try:
        # Get protein info
        original_image_path = background_remove_image_list['img_path'][protein_image_x]
        save_image_path = background_remove_image_list['df_rm_save_path'][protein_image_x]
        exposure = background_remove_image_list['value2'][protein_image_x]
        # Find nearest dark frame exposure
        df_x = parameters.find_nearest(dark_frames_list['exposure'], exposure)
        dark_frame_image_path = dark_frames_list['path'][df_x]
        # Import images
        original_image = tifffile.imread(original_image_path)
        dark_frame_image = tifffile.imread(dark_frame_image_path)

        # Subtract images
        # Get frame count
        frames = original_image.shape[0]
        frames = range(0, frames)
        # Subtract
        dark_frame_removed_image = []
        for t in frames:
            df_rm = original_image[t]
            # df_rm = df_rm - np.minimum(df_rm, dark_frame_image)
            df_rm = cv.subtract(df_rm, dark_frame_image)
            dark_frame_removed_image.append(df_rm)
        # Save image
        tifffile.imsave(save_image_path, dark_frame_removed_image, bigtiff=True, compress=tiff_compression_level, dtype=dark_frame_removed_image[0].dtype)
        time.sleep(5)
    except:
        print("Cannot complete remove")
    return original_image_path
