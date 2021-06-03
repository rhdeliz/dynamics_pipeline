import os, csv, datetime
import pandas as pd
# To prevent truncating paths
pd.options.display.max_colwidth = 10000

# Temporary files path
def processing_paths(directories_table):
    try:
        # Path
        directory_list = pd.read_csv(directories_table)

        # Input path
        input_path = directory_list.loc[directory_list['contains'] == 'input']['path']
        input_path = input_path.to_string(index=False)
        # Create if it doesn't exist
        if not os.path.exists(input_path):
            os.makedirs(input_path)

        # Processing path
        processing_path = directory_list.loc[directory_list['contains'] == 'processing']['path']
        processing_path = processing_path.to_string(index=False)
        # Create if it doesn't exist
        if not os.path.exists(processing_path):
            os.makedirs(processing_path)

        # From proprietary file to TIFF output path
        to_tiff_path = os.path.join(processing_path, "01_Convert_to_TIFF")
        # Create if it doesn't exist
        if not os.path.exists(to_tiff_path):
            os.makedirs(to_tiff_path)

        # Background remove output path
        background_remove_path = os.path.join(processing_path, "02_Remove_Backgrounds")
        # Create if it doesn't exist
        if not os.path.exists(background_remove_path):
            os.makedirs(background_remove_path)

        # Segmentation output path
        segmentation_path = os.path.join(processing_path, "03_Segment")
        # Create if it doesn't exist
        if not os.path.exists(segmentation_path):
            os.makedirs(segmentation_path)

        # Tracking output path
        tracking_path = os.path.join(processing_path, "04_Track")
        # Create if it doesn't exist
        if not os.path.exists(tracking_path):
            os.makedirs(tracking_path)

        # Final output path
        output_path = directory_list.loc[directory_list['contains'] == 'output']['path']
        output_path = output_path.to_string(index=False)
        # Create if it doesn't exist
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        # Dark frames path
        dark_frames_path = directory_list.loc[directory_list['contains'] == 'dark_frames']['path']
        dark_frames_path = dark_frames_path.to_string(index=False)

        # Flat field paths
        flat_fields_path = directory_list.loc[directory_list['contains'] == 'flat_fields']['path']
        flat_fields_path = flat_fields_path.to_string(index=False)

        # ImageJ/FIJI path
        imagej = directory_list.loc[directory_list['contains'] == 'imagej']['path']
        imagej = imagej.to_string(index=False)

        return input_path, to_tiff_path, background_remove_path, segmentation_path, tracking_path, \
               output_path, dark_frames_path, flat_fields_path, imagej
    except:
        # Get file path
        sample_path = os.path.dirname(directories_table_path)
        sample_table = 'sample ' + os.path.basename(directories_table_path)
        sample_directories_table = os.path.join(sample_path, sample_table)
        print('Input not accepted. Sample will be at ' + sample_directories_table)
        # Create table
        with open(sample_directories_table, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['contains', 'path'])
            writer.writerow(['input', 'path'])
            writer.writerow(['processing', 'path'])
            writer.writerow(['output', 'path'])
            writer.writerow(['dark_frames', 'path'])
            writer.writerow(['flat_fields', 'path'])
            writer.writerow(['ImageJ', 'path'])

# Import dark frame parameters
def dark_frame_parameters(dark_frames_table, dark_frames_path):
    try:
        # Import dark frames table
        dark_frames_list = pd.read_csv(dark_frames_table)
        # Make paths from dark-frame images
        n_df_images = len(dark_frames_list)
        n_df_images = range(0, n_df_images)
        n_df_images = list(n_df_images)
        # Combine directory with image name
        dark_frames_list = dark_frames_list.assign(
            path=lambda dataframe: dataframe['image'].map(lambda image: os.path.join(dark_frames_path, image)),
            exposure=lambda dataframe: dataframe['exposure'].map(
                lambda exposure: exposure.split(' ')[0] if exposure.split(' ')[1] == 'ms' else exposure),
            date=lambda dataframe: dataframe['image'].map(lambda image: image[:8])
        )
        # Date to python format
        dark_frames_list = dark_frames_list.assign(
            date=lambda dataframe: dataframe['date'].map(lambda date: datetime.datetime.strptime(date, '%Y%m%d'))
        )
        return dark_frames_list
    except:
        # Get file path
        sample_path = os.path.dirname(dark_frames_table)
        sample_table = 'sample ' + os.path.basename(dark_frames_table)
        sample_dark_frames_table = os.path.join(sample_path, sample_table)
        print('Input not accepted. Sample will be at ' + sample_dark_frames_table)
        # Create table
        with open(sample_dark_frames_table, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['image', 'exposure'])
            writer.writerow(['yyyymmdd img_name', 'ms'])

# Import folder list
def flat_field_parameters(flat_fields_table, flat_fields_path):
    try:
        # Import data frame list
        flat_fields_list = pd.read_csv(flat_fields_table)
        # Make paths from dark-frame images
        n_ff_images = len(flat_fields_list)
        n_ff_images = range(0, n_ff_images)
        n_ff_images = list(n_ff_images)

        # Combine directory with image name
        flat_fields_list = flat_fields_list.assign(
            path=lambda dataframe: dataframe['image'].map(lambda image: os.path.join(flat_fields_path, image)),
            date=lambda dataframe: dataframe['image'].map(lambda image: image[:8])
        )
        # Date to python format
        flat_fields_list = flat_fields_list.assign(
            date=lambda dataframe: dataframe['date'].map(lambda date: datetime.datetime.strptime(date, '%Y%m%d'))
        )

        return flat_fields_list
    except:
        # Get file path
        sample_path = os.path.dirname(flat_fields_table)
        sample_table = 'sample ' + os.path.basename(flat_fields_table)
        sample_flat_fields_table = os.path.join(sample_path, sample_table)
        print('Input not accepted. Sample will be at ' + sample_flat_fields_table)
        # Create table
        with open(sample_flat_fields_table, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['image', 'channel', 'power', 'exposure'])
            writer.writerow(['yyyymmdd img_name', 'channel', 'pct', 'ms'])

# Import images list
def image_parameters(images_table):
    try:
        # Path
        images_list = pd.read_csv(images_table)
        # Make paths from dark-frame images
        n_images = len(images_list)
        n_images = range(0, n_images)
        n_images = list(n_images)

        # Combine directory with image name
        images_list = images_list.assign(
            image=lambda dataframe: dataframe['image'].map(lambda image: os.path.splitext(image)[0]),
            date=lambda dataframe: dataframe['image'].map(lambda image: image[:8])
        )
        # Relative path
        images_list = images_list.assign(
            relative_path=images_list.apply(lambda dataframe: os.path.join(dataframe['cohort'], dataframe['image']), axis=1),
            date=lambda dataframe: dataframe['date'].map(lambda date: datetime.datetime.strptime(date, '%Y%m%d'))

        )
        return images_list, n_images
    except:
        # Get file path
        sample_path = os.path.dirname(images_table)
        sample_table = 'sample ' + os.path.basename(images_table)
        sample_images_table = os.path.join(sample_path, sample_table)
        print('Input not accepted. Sample will be at ' + sample_images_table)
        # Create table
        with open(sample_images_table, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['image', 'cohort', 'segment_with', 'ligand', 'trackmate_max_link_distance',
                             'channel protein_name', 'channel trackmate_threshold', 'channel trackmate frame gap'])
            writer.writerow(
                ['yyyymmdd img_name', 'cell_line drug parameters', 'protein_name', 'x nM search_terms', 'px',
                 'protein_name', 'signal_noise_theshold', 's'])

# Constants
def constants(constants_table):
    try:
        constants_list = pd.read_csv(constants_table)
        # TIFF Compression
        tiff_compression_level = constants_list.loc[constants_list['parameter'] == 'tiff_compression_level']['value']
        tiff_compression_level = tiff_compression_level.to_string(index=False)
        tiff_compression_level = int(tiff_compression_level)

        # Cell diameter
        cell_diameter = constants_list.loc[constants_list['parameter'] == 'cell_diameter']['value']
        cell_diameter = cell_diameter.to_string(index=False)
        cell_diameter = int(cell_diameter)

        # Puncta diameter
        puncta_diameter = constants_list.loc[constants_list['parameter'] == 'puncta_diameter']['value']
        puncta_diameter = puncta_diameter.to_string(index=False)
        puncta_diameter = int(puncta_diameter)

        return tiff_compression_level, cell_diameter, puncta_diameter
    except:
        # Get file path
        sample_path = os.path.dirname(constants_table)
        sample_table = 'sample ' + os.path.basename(constants_table)
        sample_constants_table = os.path.join(sample_path, sample_table)
        print('Input not accepted. Sample will be at ' + sample_constants_table)
        # Create table
        with open(sample_constants_table, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['parameter', 'value', 'comments'])
            writer.writerow(['tiff_compression_level', '5', 'out of 10'])
            writer.writerow(['cell_diameter', '25', 'px'])
            writer.writerow(['puncta_diameter', '5', 'px'])
