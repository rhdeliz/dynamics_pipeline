import os
import pandas as pd


def segmentation_list(images_list, file_ending, background_remove_path):
    n_images = range(0, len(images_list))
    # Add to_tiff_path
    segmentation_image_list = []
    for image_x in n_images:
        segmentation_image = images_list['segment_with'][image_x] + file_ending

        if not images_list['cohort'][image_x] == 'Calibrations':
            segmentation_image = os.path.join(background_remove_path, images_list['relative_path'][image_x],
                                              segmentation_image)
            if os.path.exists(segmentation_image):
                segmentation_image_list.append(segmentation_image)

    return segmentation_image_list






def segment_images_list(segmentation_input_path):
    # Get metadata in folder
    metadata_path = os.path.join(segmentation_input_path, 'sewlect_metadata.csv')
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


