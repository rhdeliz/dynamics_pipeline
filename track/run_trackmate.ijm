image_path = protein_path
trackmate_output_path = "protein.xml"
puncta_diameter
trackmate_threshold
trackmate_max_link_distance
trackmate_max_gap_distance = trackmate_max_link_distance * trackmate_frame_gap
trackmate_frame_gap

//Run TrackMate
run('TrackMate', "use_gui=false "
+ "save_to=["+ trackmate_output_path + "] "
+ "display_results=false "
+ "radius="+ puncta_diameter +" "
+ "threshold=" + trackmate_threshold + " "
+ "subpixel=true "
+ "median=false "
+ "channel=1 "
+ "max_distance="+ trackmate_max_link_distance +" "
+ "max_gap_distance="+ trackmate_max_gap_distance +" "
+ "max_frame_gap="+ trackmate_frame_gap +" " );
