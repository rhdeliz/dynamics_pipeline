'''
Converts ND2 file to multi-page TIFF
'''

import os, shutil, ray, time, tifffile, csv, warnings
from pims_nd2 import ND2_Reader  # Install as pims-nd2


def flatten_dictionary(dd, separator='_', prefix=''):
    return {prefix + separator + k if prefix else k: v
            for kk, vv in dd.items()
            for k, v in flatten_dictionary(vv, separator, kk).items()
            } if isinstance(dd, dict) else {prefix: dd}

@ray.remote
# images_list
def nd2_to_tiff(image_x_path, image_x_cohort, image_x_channels, to_tiff_path, tiff_compression_level):
    # Get image parameters
    try:
        image_x_name = os.path.basename(image_x_path)
        image_x_name = os.path.splitext(image_x_name)[0]
        save_path = os.path.join(to_tiff_path, image_x_cohort, image_x_name)

        if not os.path.exists(save_path):
            os.makedirs(save_path)

    except:
        print('Image parameters error')

    try:
        warnings.filterwarnings("ignore")
        img = ND2_Reader(image_x_path)
        warnings.filterwarnings("default")
    except:
        print('Cannot open ND2 image')

    try:
        # Get channel names
        n_channels = img.metadata['plane_count'] - 1

        channel_names = []
        i = 0
        while i <= n_channels:
            plane_x = 'plane_' + str(i)
            name = img.metadata[plane_x]['name']
            channel_names.append(name)
            i += 1

        channel_numbers = range(0, img.metadata['plane_count'])
        channel_numbers = list(channel_numbers)
        # Get frame count
        n_frames = img.sizes['t'] - 1
    except:
        print('Cannot get channel parameters')

    try:
        # Save metadata to text file
        file = os.path.join(save_path, 'metadata.txt')
        file = open(file, 'w')
        file.write(img.metadata_text)
        file.close()

        # Save metadata to csv
        img_metadata = flatten_dictionary(img.metadata)
        with open(os.path.join(save_path, 'metadata.csv'), 'w') as f:
            f.write('parameter, value\n')
            for key in img_metadata.keys():
                f.write('%s,%s\n' % (key, img_metadata[key]))

        # Get channel emmisions
        emissions = []
        for channel_number in channel_numbers:
            emission = 'plane_' + str(channel_number) + '_emission_nm'
            emission = img_metadata[emission]
            emission = str(int(emission))
            emissions.append(emission)

        # Get select image metadata
        search_terms = ['\r\n  N-STORM Angle: ', '\r\n  N-STORM Direction: ', '\r\n  N-STORM Focus: ', ]
        tirf_parameters = []
        for term in search_terms:
            tirf_parameter = img.metadata_text.split(term)[1]
            tirf_parameter = tirf_parameter.split('\r\n')[0]
            tirf_parameter = tirf_parameter.split(' ')[0]
            tirf_parameters.append(tirf_parameter)

        # Get illumination data
        exposures = []
        powers = []
        excitations = []
        for channel_name in channel_names:
            try:
                # Get metadata
                separator = '\r\n Name: ' + channel_name
                channel_data = img.metadata_text.split(separator)[1]
                # Get exposure
                exposure = channel_data.split('\r\n  Exposure: ')[1]
                exposure = exposure.split('\r\n')[0]
                exposures.append(exposure)
                # Get power
                power = channel_data.split('; On; ')[0]
                power = power.split(' Power:')
                power = power[len(power)-1]
                power = power.replace(' ','')
                powers.append(power)
                # Get excitation
                excitation = channel_data.split('; Power:')[0]
                excitation = excitation.split('; ExW:')
                excitation = excitation[len(excitation)-1]
                excitation = excitation.replace(' ','')
                excitations.append(excitation)
            except:
                print('Channel metadata error')
        # Get protein names
        protein_names = []
        for channel_name in channel_names:
            protein_name = image_x_channels[channel_name]
            protein_names.append(protein_name)
        # Save channel metatata
        channel_metadata = os.path.join(save_path, 'channels_metadata.csv')
        with open(channel_metadata, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['protein_name', 'channel', 'power', 'excitation', 'emmision', 'exposure', 'angle', 'direction', 'focus'])
            for c in channel_numbers:
                writer.writerow([protein_names[c], channel_names[c], powers[c], excitations[c], emissions[c], exposures[c], tirf_parameters[0], tirf_parameters[1], tirf_parameters[2]])

    except:
        print('Metadata error')

    # Split channels
    for channel_number in channel_numbers:
        try:
            frames = range(0, n_frames)
            img_channel = []
            for t in frames:
                orig = img.get_frame_2D(c=channel_number, x=0, y=0, t=t)
                img_channel.append(orig)
            # Save file
            out_name = protein_names[channel_number] + '.tif'
            out_path = os.path.join(save_path, out_name)
            tifffile.imsave(out_path, img_channel, bigtiff=True, compress=tiff_compression_level,
                            dtype=img_channel[0].dtype)
        except:
            print('Cannot split channels')

    try:
        # Move ND2 file once finished
        shutil.move(image_x_path, save_path)
        # Close image
        img.close()
        time.sleep(5)
    except:
        print('Cannot move processed nd2 file')

    return image_x_name
