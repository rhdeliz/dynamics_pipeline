import os
import numpy as np
import pandas as pd

def conversion_list(image_x, to_tiff_path, images_list):
    try:
        # Path to images and select_metadata
        img_path = os.path.join(to_tiff_path, images_list['cohort'][image_x], images_list['image_name'][image_x])
        select_metadata_path = os.path.join(img_path, 'select_metadata.csv')
        # Import table
        select_metadata = pd.read_csv(select_metadata_path)
        # Keep only channel rows
        select_metadata_filter = select_metadata['parameter1'] == 'channel'
        select_metadata = select_metadata[select_metadata_filter]
        # Reset index
        select_metadata = select_metadata.reset_index()
        # Get number of channels
        n_channels = len(select_metadata)
        n_channels = range(0, n_channels)
        n_channels = list(n_channels)

        # Get image paths
        paths = []
        for channel in n_channels:
            # Get image name
            protein_img = select_metadata['value1'][channel]
            protein_img = protein_img + '.tif'
            # Add path
            path = os.path.join(img_path, protein_img)
            paths.append(path)
        # Add paths to table
        select_metadata['img_path'] = paths
        select_metadata = select_metadata.drop(columns='index')

        # Make DFRm save names
        save_paths = []
        for channel in n_channels:
            # Get image name
            protein_img = select_metadata['value1'][channel]
            protein_img = protein_img + '_darkframe_removed.tif'
            # Add path
            path = os.path.join(img_path, protein_img)
            save_paths.append(path)
        # Add paths to table
        select_metadata['df_rm_save_path'] = save_paths

        # Make Cell MedRm save names
        save_paths = []
        for channel in n_channels:
            # Get image name
            protein_img = select_metadata['value1'][channel]
            protein_img = protein_img + '_cell_median_removed.tif'
            # Add path
            path = os.path.join(img_path, protein_img)
            save_paths.append(path)
        # Add paths to table
        select_metadata['cell_med_rm_save_path'] = save_paths

        # Make Puncta MedRm save names
        save_paths = []
        for channel in n_channels:
            # Get image name
            protein_img = select_metadata['value1'][channel]
            protein_img = protein_img + '_puncta_median_removed.tif'
            # Add path
            path = os.path.join(img_path, protein_img)
            save_paths.append(path)
        # Add paths to table
        select_metadata['puncta_med_rm_save_path'] = save_paths

    except:
        print("Cannot complete conversion_list")
    return select_metadata

# Get all TIFF images
def get_conversion_list(n_images, to_tiff_path, images_list):
    try:
        background_remove_image_list = conversion_list(0, to_tiff_path, images_list)
        if max(n_images) >= 1:
            for image_x in n_images[1:]:
                select_metadata_x = conversion_list(image_x, to_tiff_path, images_list)
                background_remove_image_list = background_remove_image_list.append(select_metadata_x)

        background_remove_image_list = background_remove_image_list.reset_index()
    except:
        print("Cannot complete get_conversion_list")
    return background_remove_image_list

# Gets nearest number in array, returns index
# For finding nearest dark-frame exposure
def find_nearest(array, value):
    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return index
