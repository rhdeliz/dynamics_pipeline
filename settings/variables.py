import os
import pandas as pd # Needs xlrd, openpyxl

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

# Import folder list
def flat_fields_parameters(parameters_table):
    # Path
    directory_list = pd.read_excel(parameters_table, 'directories')
    # Dark frames
    flat_fields_path = directory_list.loc[directory_list['contains'] == 'flat_fields']['path']
    flat_fields_path = flat_fields_path.reset_index()['path']
    flat_fields_path = flat_fields_path[0]
    # Import data frame list
    flat_fields_list = pd.read_excel(parameters_table, 'flat_fields')
    # Make paths from dark-frame images
    n_ff_images = len(flat_fields_list)
    n_ff_images = range(0, n_ff_images)
    n_ff_images = list(n_ff_images)
    # Combine directory with image name
    ff_image_path = []
    for f in n_ff_images:
        path = os.path.join(flat_fields_path, flat_fields_list['image'][f])
        ff_image_path.append(path)
    # Add loop result to list
    flat_fields_list['path'] = ff_image_path
    # Get dates
    ff_date = []
    for f in n_ff_images:
        image = flat_fields_list['image'][f]
        date = image[:8]
        ff_date.append(date)
    # Add loop result to list
    flat_fields_list['date'] = ff_date

    return flat_fields_list

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
    # TIFF Compression
    tiff_compression_level = constants_list.loc[constants_list['parameter'] == 'tiff_compression_level']['value']
    tiff_compression_level = tiff_compression_level.to_string(index=False)
    tiff_compression_level = int(tiff_compression_level)

    # Median parameters
    # Cell size for median subtraction
    cell_size = constants_list.loc[constants_list['parameter'] == 'cell_size']['value']
    cell_size = cell_size.to_string(index=False)
    cell_size = int(cell_size)
    # Puncta size for median filter and trackmate
    puncta_size = constants_list.loc[constants_list['parameter'] == 'puncta_size']['value']
    puncta_size = puncta_size.to_string(index=False)
    puncta_size = int(puncta_size)

    # Trackmate parameters
    trackmate_max_link_distance = constants_list.loc[constants_list['parameter'] == 'trackmate_max_link_distance'][
        'value']
    trackmate_max_link_distance = float(trackmate_max_link_distance)

    trackmate_max_gap_distance = constants_list.loc[constants_list['parameter'] == 'trackmate_max_gap_distance'][
        'value']
    trackmate_max_gap_distance = float(trackmate_max_gap_distance)

    trackmate_frame_gap = constants_list.loc[constants_list['parameter'] == 'trackmate_frame_gap']['value']
    trackmate_frame_gap = float(trackmate_frame_gap)

    return tiff_compression_level,cell_size, puncta_size, trackmate_max_link_distance, trackmate_max_gap_distance, trackmate_frame_gap
