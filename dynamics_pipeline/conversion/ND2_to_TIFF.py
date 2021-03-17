'''
Converts ND2 file to multi-page TIFF
'''

import os, shutil, ray, time
import conversion.channels as channels
from pims_nd2 import ND2_Reader

@ray.remote
def nd2_to_tiff(image_x, images_list, to_tiff_path, tiff_compression_level):
    # Get image parameters
    try:
        image_x_path = images_list['input_path'][image_x]
        image_x_cohort = images_list['cohort'][image_x]
        image_x_name = images_list['image_name'][image_x]
        save_path = os.path.join(to_tiff_path, image_x_cohort, image_x_name)

        if not os.path.exists(save_path):
            os.makedirs(save_path)

    except:
        print('Image parameters error')

    try:
        # Open image
        img = ND2_Reader(image_x_path)
    except:
        print('Cannot open ND2 image')

    try:
        # Get channel names
        n_channels = img.metadata['plane_count'] - 1
        channel_names = channels.get_channel_names(img, n_channels)
        channel_range = range(0, img.metadata['plane_count'])
        channel_range = list(channel_range)
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

        # Save relevant variables
        frame_rate = img.frame_rate
        frame_rate = str(frame_rate)
        pixel_size = img.calibration
        pixel_size = str(pixel_size)
        x_size = img.metadata['tile_width']
        x_size = str(x_size)
        y_size = img.metadata['tile_height']
        y_size = str(y_size)
        exposures = []
        for channel in channel_names:
            separator = '\r\n Name: ' + channel
            exposure = img.metadata_text.split(separator)[1]
            exposure = exposure.split('\r\n  Exposure: ')[1]
            exposure = exposure.split('\r\n')[0]
            exposure = exposure.split(' ')
            protein_x = images_list[channel][image_x]
            # Join
            text = 'channel,'+ protein_x + ',' + channel +',exposure,'+exposure[0]+','+exposure[1]
            exposures.append(text)
        exposures = '\n'.join(exposures)
        select_metadata =\
            'parameter1,value1,units1,parameter2,value2,units2\n'+\
            'frame_rate,' + frame_rate + ',s\n'+\
            'pixel_size,' + pixel_size + ',um\n'+ \
            'x_size,' + x_size + ',px\n' + \
            'y_size,' + y_size + ',px,\n' + \
            exposures

        # Save select metadata variables to text file
        file = os.path.join(save_path, 'select_metadata.csv')
        file = open(file, 'w')
        file.write(select_metadata)
        file.close()

    except:
        print('Metadata error')

    try:
        # Split channels
        for channel in channel_range:
            channels.save_channel_img(img, images_list, n_frames, channel, channel_names, image_x, save_path, tiff_compression_level)
    except:
        print('Cannot split channels')

    try:
        # Move ND2 file once finished
        out_name = image_x_name + '.nd2'
        move_path = os.path.join(save_path, out_name)
        shutil.move(image_x_path, move_path)
        # Close image
        img.close()
        time.sleep(3)
    except:
        print('Cannot move processed nd2 file')

    return image_x_name
