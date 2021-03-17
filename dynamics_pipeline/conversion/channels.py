import os, tifffile
import numpy as np
from PIL import Image

# Function to get channel names
def get_channel_names(img, planes):
    try:
        channels = []
        i = 0
        while i <= planes:
            plane_x = 'plane_' + str(i)
            name = img.metadata[plane_x]['name']
            channels.append(name)
            i += 1
    except:
        print('get_channel_names error')
    return channels

# Function to extract images for a channel
def get_original_img(img, frames, channel):
    try:
        t = 0
        all_orig = []
        while t <= frames:
            orig = img.get_frame_2D(c=channel, x=0, y=0, t=t)
            orig = np.array(orig)
            orig = Image.fromarray(orig)
            all_orig.append(np.array(orig))
            t += 1
    except:
        print('original error')
    return np.array(all_orig)


# Function to save splitted channel
def save_channel_img(img, images_list, n_frames, channel, channel_names, image_x, save_path, tiff_compression_level):
    # Get image
    img_channel = get_original_img(img, n_frames, channel)
    # Save file
    protein_x = images_list[channel_names[channel]][image_x]
    out_name = protein_x + '.tif'
    out_path = os.path.join(save_path, out_name)
    tifffile.imsave(out_path, img_channel, compress=tiff_compression_level)
    return out_path
