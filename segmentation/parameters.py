import os
import pandas as pd
import re


def segmentation_list(background_remove_path, segmentation_path, images_list):
    segmentation_input_path = []
    segmentation_output_path = []

    n_images = len(images_list)
    n_images = range(0, n_images)
    n_images = list(n_images)

    for image_x in n_images:
        # Input
        image_x_input_name = images_list['segment_with'][image_x] + '_cell_median_removed.tif'
        image_x_input_path = os.path.join(background_remove_path, images_list['cohort'][image_x],
                                          images_list['image_name'][image_x], image_x_input_name)
        segmentation_input_path.append(image_x_input_path)

        # Output path
        image_x_output_path = os.path.join(segmentation_path, images_list['cohort'][image_x],
                                           images_list['image_name'][image_x])
        segmentation_output_path.append(image_x_output_path)
        # Delete temporary variables
        del image_x_input_path, image_x_output_path
    # Make it a table
    segmentation_image_list = \
        {'input_path': segmentation_input_path,
         'output_path': segmentation_output_path}

    # Convert it to dataframe
    segmentation_image_list = pd.DataFrame(segmentation_image_list)
    # Delete temporary variables
    del segmentation_input_path, segmentation_output_path

    # Take out calibration images
    n_segmentation_images = images_list.image_name.str.contains('calibration', flags=re.IGNORECASE)
    n_segmentation_images = n_segmentation_images.index[n_segmentation_images == False].tolist()

    return segmentation_image_list, n_segmentation_images

def segment_images_list(segmentation_input_path):
    # Get metadata in folder
    metadata_path = os.path.join(segmentation_input_path, 'select_metadata.csv')
    select_metadata = pd.read_csv(metadata_path)
    keep_rows = select_metadata.parameter1 == 'channel'
    keep_rows = select_metadata.index[keep_rows == True].tolist()
    select_metadata = select_metadata.iloc[keep_rows]
    # List of cell images
    image_to_segment = select_metadata.value1
    image_to_segment_1 = image_to_segment + '_cell_median_removed.tif'
    image_to_segment_1 = list(image_to_segment_1)
    # List of puncta images
    image_to_segment_2 = image_to_segment + '_puncta_median_removed.tif'
    image_to_segment_2 = list(image_to_segment_2)
    image_to_segment = image_to_segment_1 + image_to_segment_2

    image_to_segment_path = []
    for img in image_to_segment:
        img_path = os.path.join(segmentation_input_path, img)
        image_to_segment_path.append(img_path)

    # ROI image
    segmentation_image = os.path.join(segmentation_input_path, 'Segmentation.tif')
    segmentation_image = [segmentation_image]*len(image_to_segment_path)

    images_to_segment = pd.DataFrame([image_to_segment_path, segmentation_image])
    images_to_segment = images_to_segment.transpose()
    images_to_segment.columns = ['image_path', 'segmentation_image_path']

    return images_to_segment


