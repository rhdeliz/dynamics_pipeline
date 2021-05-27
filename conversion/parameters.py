import os, tifffile

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
# May not need PIL?
def get_original_img(img, n_frames, channel):
    try:
        frames = range(0, n_frames)
        all_orig = []
        for t in frames:
            orig = img.get_frame_2D(c=channel, x=0, y=0, t=t)
            all_orig.append(orig)
    except:
        print('original error')
    return all_orig

# Function to save splitted channel
def save_channel_img(img, images_list, n_frames, channel, channel_names, image_x, save_path, tiff_compression_level):
    # Get image
    img_channel = get_original_img(img, n_frames, channel)
    # Save file
    protein_x = images_list[channel_names[channel]][image_x]
    out_name = protein_x + '.tif'
    out_path = os.path.join(save_path, out_name)
    tifffile.imsave(out_path, img_channel, bigtiff=True, compress=tiff_compression_level, dtype=img_channel[0].dtype)
    return out_path