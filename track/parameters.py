import os
import pandas as pd


def tracking_list(images_list, segmentation_path, input_path):
    n_images = range(0, len(images_list))
    all_channels_metadata = []
    for image_x in n_images:
        # Get paths
        relative_path = images_list['relative_path'][image_x]
        image_path = os.path.join(segmentation_path, relative_path)
        image_name = images_list['image'][image_x]
        cohort_name = images_list['cohort'][image_x]

        # Get image parameters
        ligand = images_list['ligand'][image_x]
        trackmate_max_link_distance = images_list['trackmate_max_link_distance'][image_x]

        # Get tables to merge
        select_images_list = images_list['image'] == image_name
        select_images_list = list(select_images_list)
        select_images_list = images_list.loc[select_images_list]
        # Get table paths
        channels_csv_path = os.path.join(image_path, 'channels_metadata.csv')
        area_csv_path = os.path.join(image_path, 'cell_area.csv')
        metadata_csv_path = os.path.join(image_path, 'metadata.csv')
        if os.path.exists(channels_csv_path) and os.path.exists(area_csv_path) and os.path.exists(metadata_csv_path):
            # Get tables
            metadata = pd.read_csv(metadata_csv_path)
            width = metadata.loc[metadata['parameter'] == 'width']['value']
            width = int(width)
            height = metadata.loc[metadata['parameter'] == 'height']['value']
            height = int(height)
            calibration_um = metadata.loc[metadata['parameter'] == 'calibration_um']['value']
            calibration_um = float(calibration_um)
            time_start = metadata.loc[metadata['parameter'] == 'time_start']['value']
            time_start = time_start.to_string(index=False)
            objective = metadata.loc[metadata['parameter'] == 'objective']['value']
            objective = objective.to_string(index=False)
            frame_rate = metadata.loc[metadata['parameter'] == 'frame_rate']['value']
            frame_rate = float(frame_rate)

            # Get channel metadata
            channel_metadata = pd.read_csv(channels_csv_path)
            cells = pd.read_csv(area_csv_path)
            image_channels_metadata = []
            # Get channel parameters
            for channel_x in range(0, len(channel_metadata)):
                # Channel metadata
                channel_name = channel_metadata['channel'][channel_x]
                protein_name = channel_metadata['protein_name'][channel_x]
                select_channel_metadata = channel_metadata.loc[[channel_x]]

                # Get trackmate parameters
                trackmate_threhsold = channel_name + ' trackmate_threshold'
                trackmate_threhsold = select_images_list[trackmate_threhsold]
                trackmate_threhsold = float(trackmate_threhsold)
                select_channel_metadata['trackmate_threhsold'] = trackmate_threhsold

                trackmate_frame_gap = channel_name + ' trackmate_frame_gap'
                trackmate_frame_gap = select_images_list[trackmate_frame_gap]
                trackmate_frame_gap = int(trackmate_frame_gap)
                select_channel_metadata['trackmate_frame_gap'] = trackmate_frame_gap

                trackmate_gap_link_distance = trackmate_max_link_distance
                select_channel_metadata['trackmate_gap_link_distance'] = trackmate_gap_link_distance

                # Get protein image paths
                protein_relative_paths = []
                cell_areas = []
                position_x_list = []
                position_y_list = []
                for cell_x in range(0, len(cells)):
                    cell_name = cells['cell'][cell_x]

                    protein_relative_path = os.path.join(relative_path, cell_name, protein_name)
                    protein_relative_paths.append(protein_relative_path)

                    cell_area = cells['area'][cell_x]
                    cell_areas.append(cell_area)

                    position_x = cells['position_x'][cell_x]
                    position_x_list.append(position_x)

                    position_y = cells['position_y'][cell_x]
                    position_y_list.append(position_y)


                # Expand rows and incorporate cell data
                select_channel_metadata = select_channel_metadata.loc[select_channel_metadata.index.repeat(len(cells))]
                select_channel_metadata['protein_relative_path'] = protein_relative_paths
                select_channel_metadata['area'] = cell_areas
                select_channel_metadata['position_x'] = position_x_list
                select_channel_metadata['position_y'] = position_y_list
                select_channel_metadata['ligand'] = ligand
                select_channel_metadata['cohort'] = cohort_name
                select_channel_metadata['trackmate_max_link_distance'] = trackmate_max_link_distance
                select_channel_metadata['width'] = width
                select_channel_metadata['height'] = height
                select_channel_metadata['calibration_um'] = calibration_um
                select_channel_metadata['time_start'] = time_start
                select_channel_metadata['objective'] = objective
                select_channel_metadata['frame_rate'] = frame_rate

                all_channels_metadata.append(select_channel_metadata)
                image_channels_metadata.append(select_channel_metadata)

            # Save
            image_channels_metadata = pd.concat(image_channels_metadata)
            csv_name = os.path.join(image_path, 'summmary.csv')
            image_channels_metadata.to_csv(csv_name, index = False)
    # Concatenate
    all_channels_metadata = pd.concat(all_channels_metadata)
    all_channels_metadata = all_channels_metadata.reset_index()
    all_channels_metadata= all_channels_metadata.drop(columns=['index'])
    csv_name = os.path.join(input_path, 'summary.csv')
    all_channels_metadata.to_csv(csv_name, index=False)

    return all_channels_metadata

def tracking_parameters(all_channels_metadata, file_ending, segmentation_path):
    n_trackings = range(0, len(all_channels_metadata))

    protein_paths = []
    image_paths = []
    trackmate_threhsolds = []
    trackmate_frame_gaps = []
    trackmate_max_link_distances = []
    trackmate_gap_link_distances = []

    for row in n_trackings:
        relative_path = all_channels_metadata['protein_relative_path'][row]

        protein_path = os.path.join(segmentation_path, relative_path)
        protein_paths.append(protein_path)

        image_path = protein_path + file_ending
        image_paths.append(image_path)

        trackmate_threhsold = all_channels_metadata['trackmate_threhsold'][row]
        trackmate_threhsolds.append(trackmate_threhsold)

        trackmate_frame_gap = all_channels_metadata['trackmate_frame_gap'][row]
        trackmate_frame_gaps.append(trackmate_frame_gap)

        trackmate_max_link_distance = all_channels_metadata['trackmate_max_link_distance'][row]
        trackmate_max_link_distances.append(trackmate_max_link_distance)

        trackmate_gap_link_distance = all_channels_metadata['trackmate_gap_link_distance'][row]
        trackmate_gap_link_distances.append(trackmate_gap_link_distance)

    return n_trackings, image_paths, protein_paths, trackmate_threhsolds, trackmate_frame_gaps, trackmate_max_link_distances, trackmate_gap_link_distances
