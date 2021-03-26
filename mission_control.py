import ray
import settings.parallel
import settings.variables as variables
import conversion.ND2_to_TIFF as convert_files
import background.parameters
import background.darkframe
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
to_tiff_path, background_remove_path, segmentation_path, tracking_path, output_path = variables.processing_paths(parameters_table)
median_blur_size, tiff_compression_level, median_filter_size, trackmate_blob_diameter, trackmate_max_link_distance, trackmate_max_gap_distance, trackmate_frame_gap = variables.constants(parameters_table)

# Convert nd2 to tiff
ray.init()
result_ids = []
for image_x in n_images:
    result_id = convert_files.nd2_to_tiff.remote(image_x, images_list, to_tiff_path, tiff_compression_level)
    result_ids.append(result_id)
results = settings.parallel.ids_to_vals(result_ids)
print(results)
ray.shutdown()

# Subtract dark frame
# Get list of images to have background subtracted
background_remove_image_list = background.parameters.get_conversion_list(n_images, to_tiff_path, images_list)
# Get number of protein images
n_protein_images = len(background_remove_image_list)
n_protein_images = range(0, n_protein_images)
# Run in parallel
ray.init()
result_ids = []
for protein_image_x in n_protein_images:
    result_id = background.darkframe.remove.remote(protein_image_x, background_remove_image_list, dark_frames_list, tiff_compression_level)
    result_ids.append(result_id)

results = settings.parallel.ids_to_vals(result_ids)
print(results)
ray.shutdown()

