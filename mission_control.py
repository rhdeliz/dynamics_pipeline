import ray
import settings.parallel
import settings.variables as variables

'''
conda create --name ray
conda activate ray
conda install --name ray pip
pip install ray
'''
# PIL install as Pillow
# skimage


# Parameters table
parameters_table = '/Users/u_deliz/Desktop/NewPipeline/Input/parameters.xlsx'
dark_frames_list = variables.dark_frames_parameters(parameters_table)
flat_fields_list = variables.flat_fields_parameters(parameters_table)
images_list, n_images = variables.images_parameters(parameters_table)
to_tiff_path, background_remove_path, segmentation_path, tracking_path, output_path = variables.processing_paths(
    parameters_table)
tiff_compression_level, cell_size, puncta_size, trackmate_max_link_distance, trackmate_max_gap_distance, trackmate_frame_gap = variables.constants(
    parameters_table)

# Convert nd2 to tiff
import conversion.ND2_to_TIFF as convert_files

ray.init()
result_ids = []
for image_x in n_images:
    result_id = convert_files.nd2_to_tiff.remote(image_x, images_list, to_tiff_path, tiff_compression_level)
    result_ids.append(result_id)
results = settings.parallel.ids_to_vals(result_ids)
print(results)
ray.shutdown()

# Subtract dark frame
import background.parameters
import background.darkframe

# Get list of images to have background subtracted
background_remove_image_list = background.parameters.get_conversion_list(n_images, to_tiff_path, images_list)
# Get number of protein images
n_protein_images = len(background_remove_image_list)
n_protein_images = range(0, n_protein_images)
# Run in parallel
ray.init()
result_ids = []
for protein_image_x in n_protein_images:
    result_id = background.darkframe.remove_frame.remote(protein_image_x, background_remove_image_list,
                                                         dark_frames_list, tiff_compression_level)
    result_ids.append(result_id)
results = settings.parallel.ids_to_vals(result_ids)
print(results)
ray.shutdown()

# Subtract median
import background.parameters
import background.median

# Get list of images to have background subtracted
background_remove_image_list = background.parameters.get_conversion_list(n_images, to_tiff_path, images_list)
# Get number of protein images
n_protein_images = len(background_remove_image_list)
n_protein_images = range(0, n_protein_images)
# Run
for protein_image_x in n_protein_images:
    background.median.remove_median_blur(protein_image_x, background_remove_image_list, cell_size, puncta_size,
                                         tiff_compression_level)
# Move directories of files with removed background
import os
import shutil

for image_x in n_images:
    image_x_input_path = os.path.join(to_tiff_path, images_list['cohort'][image_x], images_list['image_name'][image_x])
    image_x_output_path = os.path.join(background_remove_path, images_list['cohort'][image_x],
                                       images_list['image_name'][image_x])
    shutil.move(image_x_input_path, image_x_output_path)
    del image_x_input_path, image_x_output_path


# Segment images
import segmentation.parameters
import segmentation.operations
import background.parameters
import os


segmentation_image_list, n_segmentation_images = segmentation.parameters.segmentation_list(background_remove_path, segmentation_path, images_list)
# Create folder if it doesn't exist
if not os.path.exists(segmentation_path):
    os.mkdir(segmentation_path)

# Run in parallel
ray.init()
result_ids = []
for image_x in n_segmentation_images:
    # Get image parameters
    img_parameters = segmentation_image_list.loc[image_x]
    segmentation_input_path = img_parameters['input_path']
    # Run segmentation
    result_id = segmentation.operations.segment(segmentation_input_path, cell_size, tiff_compression_level)
    result_ids.append(result_id)
results = settings.parallel.ids_to_vals(result_ids)
print(results)
ray.shutdown()

# Make substacks
import pandas as pd
images_to_segment = []
for image_x in n_images:
    # Get image parameters
    img_parameters = segmentation_image_list.loc[image_x]
    segmentation_input_path = img_parameters['input_path']
    segmentation_input_path = os.path.dirname(segmentation_input_path)
    # Get image parameters
    temp_images_to_segment = segmentation.parameters.segment_images_list(segmentation_input_path)
    images_to_segment.append(temp_images_to_segment)
# Combine list and make it clean
images_to_segment = pd.concat(images_to_segment)
images_to_segment = images_to_segment.reset_index()
images_to_segment = images_to_segment.drop(columns=['index'])

# Run in parallel
ray.init()
result_ids = []
n_segmentations = range(0, len(images_to_segment))
n_segmentations = list(n_segmentations)
for image_x in n_segmentations:
    # Get image parameters
    img_parameters = images_to_segment.loc[image_x]
    image_x_path = img_parameters['image_path']
    segmentation_image_x_path = img_parameters['segmentation_image_path']
    # Run segmentation
    result_id = segmentation.operations.make_substacks.remote(image_x_path, segmentation_image_x_path, puncta_size, tiff_compression_level)
    result_ids.append(result_id)
results = settings.parallel.ids_to_vals(result_ids)
print(results)
ray.shutdown()

# Move directories of files with segmentation
import os
import shutil

for image_x in n_images:
    image_x_input_path = os.path.join(to_tiff_path, images_list['cohort'][image_x], images_list['image_name'][image_x])
    image_x_output_path = os.path.join(background_remove_path, images_list['cohort'][image_x],
                                       images_list['image_name'][image_x])
    shutil.move(image_x_input_path, image_x_output_path)
    del image_x_input_path, image_x_output_path
