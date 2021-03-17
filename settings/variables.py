import os
import pandas as pd
# Needs xlrd, openpyxl

# Import folder list
def dark_frames_parameters(parameters_table):
    # Path
    directory_list = pd.read_excel(parameters_table, 'directories')
    # Dark frames
    dark_frames_path = directory_list.loc[directory_list['contains'] == 'dark_frames']['path']
    dark_frames_path = dark_frames_path[0]
    # Import data frame list
    dark_frames_list = pd.read_excel(parameters_table, 'dark_frames')
    # Make paths from dark-frame images
    n_df_images = len(dark_frames_list)
    n_df_images = range(0, n_df_images)
    n_df_images = list(n_df_images)
    # Combine directory with image name
    df_image_path = []
    for f in n_df_images:
        path = os.path.join(dark_frames_path, dark_frames_list['image'][f])
        df_image_path.append(path)
    # Add loop result to list
    dark_frames_list['path'] = df_image_path

    return dark_frames_list

# Import images list
def images_parameters(parameters_table):
    # Path
    directory_list = pd.read_excel(parameters_table, 'directories')
    images_list = pd.read_excel(parameters_table, 'images')
    images_path = directory_list.loc[directory_list['contains'] == 'input']['path']
    images_path = images_path.to_string(index=False)
    # Make paths from dark-frame images
    n_images = len(images_list)
    n_images = range(0, n_images)
    n_images = list(n_images)
    # Combine directory with image name
    image_path = []
    for f in n_images:
        path = os.path.join(images_path, images_list['image_name'][f])
        image_path.append(path)
    # Add loop result to list
    images_list['input_path'] = image_path
    # Take out file format in image name
    base_name = []
    for f in n_images:
        name = os.path.splitext(images_list['image_name'][f])[0]
        base_name.append(name)
    # Add loop result to list
    images_list['image_name'] = base_name

    return images_list, n_images

# Temporary files path
def processing_paths(parameters_table):
    # Path
    directory_list = pd.read_excel(parameters_table, 'directories')

    processing_path = directory_list.loc[directory_list['contains'] == 'processing']['path']
    processing_path = processing_path.to_string(index=False)
    # Paths
    to_tiff_path = os.path.join(processing_path, "01_Convert_to_TIFF")
    if not os.path.exists(to_tiff_path):
        os.makedirs(to_tiff_path)

    background_remove_path = os.path.join(processing_path, "02_Remove_Backgrounds")
    if not os.path.exists(background_remove_path):
        os.makedirs(background_remove_path)

    segmentation_path = os.path.join(processing_path, "03_Segment")
    if not os.path.exists(segmentation_path):
        os.makedirs(segmentation_path)

    tracking_path = os.path.join(processing_path, "04_Track")
    if not os.path.exists(tracking_path):
        os.makedirs(tracking_path)

    # Output path
    output_path = directory_list.loc[directory_list['contains'] == 'output']['path']
    output_path = output_path.to_string(index=False)

    return to_tiff_path, background_remove_path, segmentation_path, tracking_path, output_path

# Constants
def constants(parameters_table):
    constants_list = pd.read_excel(parameters_table, 'constants')
    # Median blur size to subtract from dark-frame removed image
    median_blur_size = constants_list.loc[constants_list['parameter'] == 'median_pixel_size']['value']
    median_blur_size = median_blur_size.to_string(index=False)
    # TIFF Compression
    tiff_compression_level = constants_list.loc[constants_list['parameter'] == 'tiff_compression_level']['value']
    tiff_compression_level = tiff_compression_level.to_string(index=False)
    tiff_compression_level = int(tiff_compression_level)
    # Real median or approximate
    real_median = constants_list.loc[constants_list['parameter'] == 'real_median']['value']
    real_median = real_median.to_string(index=False)
    real_median = real_median.lower()
    real_median = real_median == 'yes'
    # Trackmate parameters
    trackmate_blob_diameter = constants_list.loc[constants_list['parameter'] == 'trackmate_blob_diameter']['value']
    trackmate_blob_diameter = float(trackmate_blob_diameter)

    trackmate_max_link_distance = constants_list.loc[constants_list['parameter'] == 'trackmate_max_link_distance'][
        'value']
    trackmate_max_link_distance = float(trackmate_max_link_distance)

    trackmate_max_gap_distance = constants_list.loc[constants_list['parameter'] == 'trackmate_max_gap_distance'][
        'value']
    trackmate_max_gap_distance = float(trackmate_max_gap_distance)

    trackmate_frame_gap = constants_list.loc[constants_list['parameter'] == 'trackmate_frame_gap']['value']
    trackmate_frame_gap = float(trackmate_frame_gap)

    return median_blur_size, tiff_compression_level, real_median, trackmate_blob_diameter, trackmate_max_link_distance, trackmate_max_gap_distance, trackmate_frame_gap
