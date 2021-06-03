import os
import pandas as pd

def get_images_to_track(segmentation_path, images_list, image_x):
    try:
        # Get image path
        substacked_image = os.path.join(segmentation_path, images_list['cohort'][image_x],
                                        images_list['image_name'][image_x])
        # Get proteins in image
        metadata_table = os.path.join(substacked_image, 'select_metadata.csv')
        metadata_table = pd.read_csv(metadata_table)
        keep_rows = metadata_table.parameter1 == 'channel'
        keep_rows = keep_rows.index[keep_rows == True].tolist()
        metadata_table = metadata_table.iloc[keep_rows]
        substacked_proteins = metadata_table['value1']
        substacked_proteins = list(substacked_proteins)

        # Get max cell number
        cell_folders = os.listdir(substacked_image)
        cell_folders = [x for x in cell_folders if 'Cell_' in x]
        # Make path of cell protein image
        cell_paths = []
        for cell in cell_folders:
            cell_path = os.path.join(substacked_image, cell)
            cell_paths.append(cell_path)

        cell_paths = cell_paths * len(substacked_proteins)
        proteins = substacked_proteins * len(cell_folders)

        cell_paths.sort()
        image_paths = []
        for x in range(0, len(cell_paths)):
            image_name = proteins[x] + '_puncta_median_removed.tif'
            image_path = os.path.join(cell_paths[x], image_name)
            image_paths.append(image_path)
        substacked_image = os.path.join(segmentation_path, images_list['cohort'][image_x],
                                        images_list['image_name'][image_x])
        # Get proteins in image
        metadata_table = os.path.join(substacked_image, 'select_metadata.csv')
        metadata_table = pd.read_csv(metadata_table)
        keep_rows = metadata_table.parameter1 == 'channel'
        keep_rows = keep_rows.index[keep_rows == True].tolist()
        metadata_table = metadata_table.iloc[keep_rows]
        substacked_proteins = metadata_table['value1']
        substacked_proteins = list(substacked_proteins)
        # Get protein channel name
        substacked_channels = metadata_table['units1']
        substacked_channels = list(substacked_channels)

        # Get max cell number
        cell_folders = os.listdir(substacked_image)
        cell_folders = [x for x in cell_folders if 'Cell_' in x]
        # Make path of cell protein image
        cell_paths = []
        for cell in cell_folders:
            cell_path = os.path.join(substacked_image, cell)
            cell_paths.append(cell_path)

        cell_paths = cell_paths * len(substacked_proteins)
        proteins = substacked_proteins * len(cell_folders)
        channels = substacked_channels * len(cell_folders)

        cell_paths.sort()
        image_paths = []
        xml_paths = []
        for x in range(0, len(cell_paths)):
            image_name = proteins[x] + '_puncta_median_removed.tif'
            image_path = os.path.join(cell_paths[x], image_name)
            image_paths.append(image_path)

            xml_name = proteins[x] + '.xml'
            xml_path = os.path.join(cell_paths[x], xml_name)
            xml_paths.append(xml_path)



        return cell_paths, image_paths, xml_paths, proteins, channels
    except:
        print("Cannot get list of cells for tracking")
        cell_paths = []
        image_paths = []
        proteins = []
        channels = []
        xml_paths = []
        return cell_paths, image_paths, xml_paths, proteins, channels

def tracking_list(segmentation_path, images_list):
    # Get number of images
    n_images = range(len(images_list))
    n_images = list(n_images)

    # Get parameters
    cell_paths = []
    image_paths = []
    xml_paths = []
    proteins = []
    channels = []
    for image_x in n_images:
        try:
            cell_path, image_path, xml_path, protein, channel = get_images_to_track(segmentation_path, images_list, image_x)
            cell_paths = cell_paths + cell_path
            image_paths = image_paths + image_path
            xml_paths = xml_paths + xml_path
            proteins = proteins + protein
            channels = channels + channel
        except:
            print('')

    # Tracking list

    images_to_track = pd.DataFrame([protein_paths, proteins, channels])
    images_to_track = images_to_track.transpose()
    images_to_track.columns = ['protein_image_path', 'protein', 'channel']

    return images_to_track

def imagej_path(parameters_table):
    # Path
    directory_list = pd.read_excel(parameters_table, 'directories')
    # Dark frames
    flat_fields_path = directory_list.loc[directory_list['contains'] == 'ImageJ']['path']
    flat_fields_path = flat_fields_path.reset_index()['path']
    flat_fields_path = flat_fields_path[0]

    return flat_fields_path

def tracking_parameters():


