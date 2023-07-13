// Macro for converting avi-videos to bitmap images, obtaining the desired command file settings and running CellTraxx within this ImageJ macro
// (C) Børge Holme, SINTEF, 2023-05-26

// Fresh start
close("*");
run("Clear Results");

// DECLARING GLOBAL VARIABLES THAT ARE ALSO USED BY FUNCTIONS Get_Dialog_Settings() AND Write_Settings_To_File()

var wound_healing_mode,					version;
var perform_flat_field_correction,		track_smoothing_iterations;	
var perform_image_shift_correction,		dat_drive_letter;
var perform_interactive_tuning,			first_part_of_results_folder_name;	
var pixel_size, 						top____crop_margin;
var gaussian_filter_radius, 			bottom_crop_margin;
var smallest_cell_diameter, 			left___crop_margin;	
var largest__cell_diameter, 			right__crop_margin;
var time_between_images, 				first_image;
var highest_cell_velocity, 				last_image;
var shortest_cell_track, 				increment;
var make_identified_cell_videos, 		tracking__dot_diameter;	
var make_matched_cell_videos,			valid_track_image_contrast;													
var make_valid_track_videos,	        scalebar_color;																										
var longest_image_shift,				cutting__cell_diameter;
var segmentation_limit,					nof_bins_in_histogram,			tuning_image_code;	
var dat_dir, sys_dir, tun_dir, dat_drive_letter;

title = "CellTraxx";
dat_dir = "C:/celltraxx_data/"						// Can be changed by user in dialog window
sys_dir = "C:/celltraxx_system/"
tun_dir = "C:/celltraxx_system/tuning_temp/"
File.makeDirectory(tun_dir); 						// Just in case the user has cleaned up and deleted this temporary directory
if(  ! File.exists(tun_dir) ) exit("Unable to create the temporary tuning directory "+tun_dir);
File.setDefaultDir(sys_dir); 
print("\n");
print("========================");
print("       Welcome to CellTraxx!");
print("========================");
print("\nWill read defaults from folder: ", sys_dir); 

setBatchMode(true); 

// DEFINING BUTTON POSITIONS AND CALCULATING WIDTHS 
font_size   =  30; 
font_vertical_shift = 0.5*font_size + 4;
inset_height= 250; 
col_1__left = 301;		col_1_right = 355;
col_2__left = 414; 		col_2_right = 437;		
col_3__left = 507; 		col_3_right = 588;
col_4__left = 619; 		col_4_right = 705;
col_5__left = 710; 		col_5_right = 800;	
updat__left = 444;		updat_right = 475; 
row_1_top   =  23;		row_1_bot   =  71;		row_1_mid =  47; 
row_2_top   =  77;		row_2_bot   = 125;		row_2_mid = 101;
row_3_top   = 131;		row_3_bot   = 179;		row_3_mid = 155;
row_4_top   = 185;		row_4_bot   = 233;		row_4_mid = 209;
tic_1_top   =  81;		tic_1_bot   = 107;
tic_2_top   = 117;		tic_2_bot   = 143;
tic_3_top   = 153;		tic_3_bot   = 179;
col_2_width = col_2_right-col_2__left+1; 
updat_width = updat_right-updat__left+1; 

// READING DEFAULT SETTINGS, LETTING USER CHANGE DESIRED VALUES, WRITING SETTINGS BACK TO FILE TO BE USED BY CELLTRAXX.EXE 
filestring = File.openAsString(sys_dir+"celltraxx_defaults.txt");
textlines = split(filestring, "\n");
i = 0;
entries=split(textlines[i++]," \t");		if( entries[3]=="yes" || entries[3]=="Yes" || entries[3]=="YES" ) wound_healing_mode             = true; else wound_healing_mode             = false; 
entries=split(textlines[i++]," \t");		if( entries[4]=="yes" || entries[4]=="Yes" || entries[4]=="YES" ) perform_flat_field_correction  = true; else perform_flat_field_correction  = false; 
entries=split(textlines[i++]," \t");		if( entries[4]=="yes" || entries[4]=="Yes" || entries[4]=="YES" ) perform_image_shift_correction = true; else perform_image_shift_correction = false; 
entries=split(textlines[i++]," \t");		if( entries[3]=="yes" || entries[3]=="Yes" || entries[3]=="YES" ) perform_interactive_tuning     = true; else perform_interactive_tuning     = false; 
entries=split(textlines[i++]," \t");		version                		= parseFloat(entries[1]);
entries=split(textlines[i++]," \t");		track_smoothing_iterations  = parseInt(  entries[3]); 
entries=split(textlines[i++]," \t");		dat_drive_letter    				   = entries[4]; 
entries=split(textlines[i++]," \t");		results_folder_name 				   = entries[3]; 
entries=split(textlines[i++]," \t");		first_part_of_results_folder_name 	   = entries[5]; 	
entries=split(textlines[i++]," \t");		pixel_size             		= parseFloat(entries[3]); 
entries=split(textlines[i++]," \t");		gaussian_filter_radius 		= parseFloat(entries[4]); 
entries=split(textlines[i++]," \t");		smallest_cell_diameter 		= parseFloat(entries[4]); 
entries=split(textlines[i++]," \t");		largest__cell_diameter 		= parseFloat(entries[4]); 
entries=split(textlines[i++]," \t");		cutting__cell_diameter	    = parseFloat(entries[4]); 
entries=split(textlines[i++]," \t");		top____crop_margin     		= parseInt(  entries[4]); 
entries=split(textlines[i++]," \t");		bottom_crop_margin     		= parseInt(  entries[4]); 
entries=split(textlines[i++]," \t");		left___crop_margin     		= parseInt(  entries[4]); 
entries=split(textlines[i++]," \t");		right__crop_margin     		= parseInt(  entries[4]); 
entries=split(textlines[i++]," \t");		time_between_images    		= parseFloat(entries[4]); 
entries=split(textlines[i++]," \t");		highest_cell_velocity  		= parseFloat(entries[4]); 
entries=split(textlines[i++]," \t");		shortest_cell_track    		= parseFloat(entries[4]); 
entries=split(textlines[i++]," \t");		first_image 				= parseInt(  entries[3]); 
entries=split(textlines[i++]," \t");		last_image  				= parseInt(  entries[3]); 
entries=split(textlines[i++]," \t");		increment   				= parseInt(  entries[3]); 
entries=split(textlines[i++]," \t");		if( entries[4]=="yes" || entries[4]=="Yes" || entries[4]=="YES" ) make_identified_cell_videos    = true; else make_identified_cell_videos    = false; 
entries=split(textlines[i++]," \t");		if( entries[4]=="yes" || entries[4]=="Yes" || entries[4]=="YES" ) make_matched_cell_videos       = true; else make_matched_cell_videos       = false; 
entries=split(textlines[i++]," \t");		if( entries[4]=="yes" || entries[4]=="Yes" || entries[4]=="YES" ) make_valid_track_videos        = true; else make_valid_track_videos        = false; 
entries=split(textlines[i++]," \t");		tracking__dot_diameter 		= parseFloat(entries[4]); 
entries=split(textlines[i++]," \t");		valid_track_image_contrast  = parseFloat(entries[4]); 
entries=split(textlines[i++]," \t");		scalebar_color = entries[2];  			
entries=split(textlines[i++]," \t");		if( entries[3]=="yes" || entries[3]=="Yes" || entries[3]=="YES" ) draw_cell_outline	             = true; else draw_cell_outline              = false; 
entries=split(textlines[i++]," \t");		if( entries[4]=="yes" || entries[4]=="Yes" || entries[4]=="YES" ) draw_cell_track_line           = true; else draw_cell_track_line           = false; 
entries=split(textlines[i++]," \t");		if( entries[4]=="yes" || entries[4]=="Yes" || entries[4]=="YES" ) write_mirror_margin_images     = true; else write_mirror_margin_images     = false; 
entries=split(textlines[i++]," \t");		if( entries[3]=="yes" || entries[3]=="Yes" || entries[3]=="YES" ) write_shifted_images           = true; else write_shifted_images           = false; 
entries=split(textlines[i++]," \t");		if( entries[4]=="yes" || entries[4]=="Yes" || entries[4]=="YES" ) write_gaussian_smoothed_images = true; else write_gaussian_smoothed_images = false; 
entries=split(textlines[i++]," \t");		if( entries[4]=="yes" || entries[4]=="Yes" || entries[4]=="YES" ) write_segmented_cell_images    = true; else write_segmented_cell_images    = false; 
entries=split(textlines[i++]," \t");		if( entries[4]=="yes" || entries[4]=="Yes" || entries[4]=="YES" ) write_cut_cell_images          = true; else write_cut_cell_images          = false; 
entries=split(textlines[i++]," \t");		if( entries[4]=="yes" || entries[4]=="Yes" || entries[4]=="YES" ) write_identified_cell_images   = true; else write_identified_cell_images   = false; 
entries=split(textlines[i++]," \t");		if( entries[4]=="yes" || entries[4]=="Yes" || entries[4]=="YES" ) write_matched_cell_images      = true; else write_matched_cell_images      = false; 
entries=split(textlines[i++]," \t");		if( entries[4]=="yes" || entries[4]=="Yes" || entries[4]=="YES" ) write_valid_track_images       = true; else write_valid_track_images       = false; 
entries=split(textlines[i++]," \t");		if( entries[5]=="yes" || entries[5]=="Yes" || entries[5]=="YES" ) keep_cells_from_previous_image = true; else keep_cells_from_previous_image = false; 
entries=split(textlines[i++]," \t");		segmentation_limit     		= parseFloat(entries[3]); 
entries=split(textlines[i++]," \t");		nof_bins_in_histogram  		= parseInt(  entries[4]); 
entries=split(textlines[i++]," \t");		tuning_image_code			= parseInt(  entries[3]); 

if( perform_flat_field_correction  ) ffc_string = "_FFC"; else  ffc_string = "";			// Addition to bitmap file names if flat field correction is on, allowing both types of images to exist 
if( scalebar_color == "No_scalebar" ) scalebar_color = ""; 									// Coding no color with "No_scalebar" in settings file to avoid read error 
if( first_part_of_results_folder_name == "#" ) first_part_of_results_folder_name = ""; 		// Coding no prefix with # in settings file to avoid read error 
default_part_of_results_folder_name = "GFR="+gaussian_filter_radius+"um_d="+smallest_cell_diameter+"-"+largest__cell_diameter+"um_cut="+cutting__cell_diameter+"um"; 

perform_interactive_tuning = true; 		// Setting tuning mode as default. 

//Get_Dialog_Settings(); 					// Displaying the input settings window on screen, accepting changes by user in settings and waiting until [OK] or [Cancel] is pressed

Dialog.create("CellTraxx Settings");		
Dialog.addCheckbox("Wound healing mode",                                          wound_healing_mode); 	    Dialog.addToSameRow();	  	Dialog.addNumber("CellTraxx version",                                                        version, 1, 3, " ");
Dialog.addCheckbox("Perform flat field correction",                    perform_flat_field_correction); 		Dialog.addToSameRow();		Dialog.addNumber("Number of track smoothing iterations", track_smoothing_iterations, 0, 3, "(0 = no smoothing)"); 			
Dialog.addCheckbox("Perform image shift correction",                  perform_image_shift_correction); 		Dialog.addToSameRow();	  	Dialog.addString("Disk drive letter for data folder X:/celltraxx_data/", dat_drive_letter, 1);	
Dialog.addCheckbox("Perform interactive tuning",                          perform_interactive_tuning);  	Dialog.addToSameRow();	  	Dialog.addString("First part of results folder name", first_part_of_results_folder_name, 17);
Dialog.addMessage( "Segmentation", 18); 											
Dialog.addNumber(  "Pixel size",                         		              pixel_size, 3, 8, "um");		Dialog.addToSameRow();		Dialog.addNumber("   Top crop margin",                   top____crop_margin, 0, 8, "pixels");	
Dialog.addNumber(  "Gaussian filter radius",                      gaussian_filter_radius, 1, 8, "um"); 	 	Dialog.addToSameRow();		Dialog.addNumber("Bottom crop margin",                   bottom_crop_margin, 0, 8, "pixels");	
Dialog.addNumber(  "Smallest cell diameter",                      smallest_cell_diameter, 1, 8, "um");		Dialog.addToSameRow();		Dialog.addNumber("  Left crop margin",                   left___crop_margin, 0, 8, "pixels");			
Dialog.addNumber(  "Largest cell diameter",                       largest__cell_diameter, 1, 8, "um");		Dialog.addToSameRow();		Dialog.addNumber(" Right crop margin",                   right__crop_margin, 0, 8, "pixels");	
Dialog.addNumber(  "Cutting cell diameter",                       cutting__cell_diameter, 1, 8, "um");			
Dialog.addMessage( "Tracking", 18); 
Dialog.addNumber(  "Time between images",                       time_between_images, 1, 8, "minutes");		Dialog.addToSameRow();		Dialog.addNumber("Analyse from image number",             first_image, 0, 8, "(0 = first)" );	
Dialog.addNumber(  "Highest cell velocity",                   highest_cell_velocity, 1, 8, "um/min" );		Dialog.addToSameRow();		Dialog.addNumber("Analyse to image number",                last_image, 0, 8, "(999 = last)");
Dialog.addNumber(  "Shortest cell track",                       shortest_cell_track, 0, 8, "images" );		Dialog.addToSameRow();		Dialog.addNumber("Image number increment",                  increment, 0, 8, "image(s)"    );	
Dialog.addMessage( "Output", 18); 
Dialog.addCheckbox("Generate videos with identified cells",              make_identified_cell_videos); 		Dialog.addToSameRow();		Dialog.addNumber("Tracking  dot diameter",            tracking__dot_diameter, 1, 8, "um"   );		
Dialog.addCheckbox("Generate videos with matched cells",                    make_matched_cell_videos); 		Dialog.addToSameRow();		Dialog.addNumber("Valid track image contrast",    valid_track_image_contrast, 1, 8, "times");	
Dialog.addCheckbox("Generate videos with valid tracks",                      make_valid_track_videos); 		Dialog.addToSameRow(); 		Dialog.addString("Scale bar color", scalebar_color, 12);					

Dialog.show();	

wound_healing_mode 					= Dialog.getCheckbox();				version 							= Dialog.getNumber();
perform_flat_field_correction  		= Dialog.getCheckbox();				track_smoothing_iterations			= Dialog.getNumber();	
perform_image_shift_correction		= Dialog.getCheckbox();				dat_drive_letter 					= Dialog.getString();
perform_interactive_tuning			= Dialog.getCheckbox();				first_part_of_results_folder_name	= Dialog.getString();	
// Segmentation
pixel_size 							= Dialog.getNumber();				top____crop_margin 					= Dialog.getNumber();
gaussian_filter_radius 				= Dialog.getNumber();				bottom_crop_margin 					= Dialog.getNumber();
smallest_cell_diameter 				= Dialog.getNumber();				left___crop_margin 					= Dialog.getNumber();	
largest__cell_diameter 				= Dialog.getNumber();				right__crop_margin 					= Dialog.getNumber();
cutting__cell_diameter 				= Dialog.getNumber();
// Tracking			
time_between_images 				= Dialog.getNumber();				first_image			 				= Dialog.getNumber();
highest_cell_velocity 				= Dialog.getNumber();				last_image  						= Dialog.getNumber();
shortest_cell_track 				= Dialog.getNumber();				increment 							= Dialog.getNumber();
// Output
make_identified_cell_videos 		= Dialog.getCheckbox();				tracking__dot_diameter 				= Dialog.getNumber();	
make_matched_cell_videos			= Dialog.getCheckbox();				valid_track_image_contrast		 	= Dialog.getNumber();													
make_valid_track_videos	            = Dialog.getCheckbox();				scalebar_color 						= Dialog.getString();																										

// Updating variables in case checkboxes were changed by user 
if( perform_image_shift_correction    ) longest_image_shift = 30;   else longest_image_shift = 0;		// Just setting a quite large default value since the maximum image shift no longer reduces the processing speed but only sets the size of the search matrix 
if( perform_flat_field_correction     ) ffc_string = "_FFC";        else ffc_string = "";				// Addition to image names   
//	if( tracking__dot_diameter == 0 ) draw_cell_center_of_mass = false; else draw_cell_center_of_mass = true; 
if( scalebar_color == "" ) scalebar_color = "No_scalebar"; 

// Updating name of results folder based on current settings and generating full name
default_part_of_results_folder_name = "GFR="+gaussian_filter_radius+"um_d="+smallest_cell_diameter+"-"+largest__cell_diameter+"um_cut="+cutting__cell_diameter+"um"; 
if( first_part_of_results_folder_name == "" ) { results_folder_name = default_part_of_results_folder_name+ffc_string;    first_part_of_results_folder_name = "#"; }
else 							   			    results_folder_name = first_part_of_results_folder_name+"_"+default_part_of_results_folder_name+ffc_string;

char_code = charCodeAt(dat_drive_letter, 0); 						// Taking the first character of the given disk drive string, in case someone would write C:
if( char_code > 90 ) char_code -= 32; 								// Converting lower case to upper case
if( char_code < 65 || char_code > 90 ) { char_code = 67; print("Invalid disk drive letter '"+dat_drive_letter+"'. Using 'C' as default."); }	 
dat_drive_letter = fromCharCode(char_code); 
dat_dir = dat_drive_letter+":/celltraxx_data/";
if( increment < 1 ) { increment = 1; print("Invalid 'Image number increment' set to 1."); }		// Just making sure the images will be incremented even if accidentally set < 1 by user

Write_Settings_To_File(); 				// Writing settings back to file to store any changes 

// CONVERTING ALL .AVI VIDEOS IN CELLTRAXX FOLDER TO SERIES OF BITMAPS IF THIS HAS NOT BEEN DONE EARLIER. GENERATING TEXT FILE WITH NAMES AND NUMBERS OF BITMAPS TO PROCESS. 
File.saveString("Bitmap names (without XXXX.bmp)		First image number	Last image number\n", sys_dir+"celltraxx_bitmaps2process.txt");		// Putting header line in file 
list = getFileList(dat_dir);

avi_count = 0;								// To allow checking below if any avi files have been found, otherwise exiting 
if( first_image < 0 ) first_image = 0; 		// Just making sure the first image to process is a non-negative number
for (i = 0; i < list.length; i++){
	j = indexOf(list[i], ".");
 	k = list[i].length; 
	if( j>0 && j<k ) {
		extension = substring(list[i], j, k);
		if( extension==".avi" || extension==".AVI" ) {
			n = first_image;
			avi_count++;
			current_last = last_image; 			// Working on a temporary copy of last image number in case this avi has fewer images 
			bmp_name = substring(list[i], 0, j)+ffc_string+"_"+IJ.pad(n,4)+".bmp";
			if( File.exists(dat_dir+bmp_name) ) {
				print("The bitmap image", bmp_name, "already exists.\nSkipping converstion of video file ", list[i]);
				do {							// Still need to generate line for this avi in celltraxx_bitmaps2process.txt, so checking images with given increment until no file exists. 
					n += increment; 
					bmp_name = substring(list[i], 0, j)+ffc_string+"_"+IJ.pad(n,4)+".bmp";
				} while( File.exists(dat_dir+bmp_name) );
				n -= increment; 
				if( current_last > n ) current_last = n; 						// Making sure that last image to process actually exists 
 				print("First Image : ", first_image, "    Increment : ", increment, "   Last Image :", current_last); 
				File.append(substring(list[i], 0, j)+ffc_string+"_"+"				"+first_image+"			"+current_last, sys_dir+"celltraxx_bitmaps2process.txt");		
			} else {
				print("Converting video file ", list[i], "to bitmaps for processing.\n"); 
				setBatchMode(true);				// Needed to not see the video image when loaded 
				run("AVI...", "open="+dat_dir+list[i]+" use convert");
				Video = getTitle();
				imagename = getTitle; 
				k = indexOf(imagename, "."); 
				imagename = substring(imagename, 0, k)+ffc_string+"_";
				getDimensions(width, height, channels, slices, frames);
				if( current_last > slices-1 ) current_last = slices-1;		 	// Making sure that last image to process actually exists 	
				File.append(imagename+"				"+first_image+"			"+current_last, sys_dir+"celltraxx_bitmaps2process.txt");		// Generating list of file names to copy into CellTrax command file and to read back by bmp2avi.ijm. Reading one fewer image since the two last in .avi are the same.

				// ROUTINE FOR PERFORMING A TYPE OF FLATFIELD CORRECTION, ASSUMING ALL PARTS OF THE BACKGROUND ARE VISIBLE MOST OF THE TIME WHEN CELLS MOVE AROUND
	   			if( perform_flat_field_correction ) {
	   				print("Performing flat field correction on video " + substring(imagename, 0, k) + ".avi"); 
	   				print("This might take several seconds..."); 
	   				run("Z Project...", "projection=Median");
	   				run("Gaussian Blur...", "sigma=2");
	   				saveAs("bmp", dat_dir + imagename + "_Background");		// Saving smoothed median image into input image folder for checking 
	   				Median = getTitle();
					imageCalculator("Divide create 32-bit stack", Video, Median);
					setMinAndMax(0.2, 1.8);			// Since the 32-bit float image has values centred on 1.0, setting the new Min-Max limits symmetrically around 1.0
					run("8-bit");
					print("Finished flat field correction of video " + substring(imagename, 0, k) + ".avi\n"); 
				} else   print("No flat field correction for video " + substring(imagename, 0, k) + ".avi"); 

//				run("Image Sequence... ", "format=BMP start=0 digits=4 name="+imagename+" save=dat_dir");		// Works in ImageJ 1.52t
				run("Image Sequence... ", "format=BMP start=0 digits=4 name="+imagename+" save="+dat_dir);		// Works in ImageJ 1.53t    Apparently a shift in policy of interpretation: https://forum.image.sc/t/following-update-macro-loading-images-as-sequence-no-longer-works/56996	&	https://imagej.nih.gov/ij/notes.html
				close();
				setBatchMode(false); 
			}
			if( avi_count == 1 ) {		// It is the first .avi file so copying the first image into tuning image to display the right image from the start
				bmp_name = substring(list[i], 0, j)+ffc_string+"_"+IJ.pad(first_image,4)+".bmp";
				print("Generating tuning image from bitmap ", bmp_name); 
				setBatchMode(true);				// Needed to not see the image when loaded 
				open(dat_dir+bmp_name); 
				Image.copy; 
				getDimensions(width, height, channels, slices, frames);
				newImage("Tuning_image", "RGB", width, height+inset_height, 1);
				Image.paste(0,inset_height);
				open(sys_dir+"celltraxx_tuning_image_inset.bmp"); 
				Image.copy; 
				selectImage("Tuning_image");
				Image.paste(0,0);
				save(tun_dir+"celltraxx_tuning_image.bmp"); 
				close("*"); 
				setBatchMode(false); 		
//				selectImage("Tuning_image");
			}
		}	
	} 
}
if( avi_count == 0 ) {
	print("No .avi files in data folder "+dat_dir);
	showMessage("Error message", "No .avi files in the data folder  "+dat_dir+"\nPlease copy the .avi videos to analyse into that folder.");
	exit;
}

File.append("stop", sys_dir+"celltraxx_bitmaps2process.txt");	
//File.close(); 

if( perform_interactive_tuning ) {
	print("Starting interactive tuning...");
	setBatchMode(false); 										// Needed to display tuning image on screen
	open(tun_dir+"celltraxx_tuning_image.bmp");
	rename("CellTraxx Tuning Image");

	x = y = 100;
	flags = 0;
	original_gaussian_filter_radius = gaussian_filter_radius;	// Storing original values in case [Reset] is clicked
	original_smallest_cell_diameter = smallest_cell_diameter;
	original_largest__cell_diameter = largest__cell_diameter;
	original_cutting__cell_diameter = cutting__cell_diameter;
	original_top____crop_margin = top____crop_margin;
	original_bottom_crop_margin = bottom_crop_margin;
	original_left___crop_margin = left___crop_margin;
	original_right__crop_margin = right__crop_margin;
	setKeyDown("none"); 
	setFont("SansSerif", font_size,"bold");
	setTool("rectangle");										// Makes sure the selected tool button is rectangle, in case some other button was selected when the macro was started 
	setColor("magenta");
	setLineWidth(1);

	// Running celltraxx.exe once with the default settings to show the effect they have 
	exec("cmd /c start "+sys_dir+"celltraxx");    					// "start" will open DOS window for output during execution 
	filestring = File.openAsString(sys_dir+"celltraxx_error_messages.txt");	// Checking if celltraxx.exe left a warning or error message, if so exiting and writing the error text to a message box
	if(      lengthOf(filestring) > 11 ) { close("CellTraxx Tuning Image"); exit(filestring); }  	// Exiting macro after closing tuning window and writing the error message from celltraxx.exe to message box
	else if( lengthOf(filestring) == 0 ) { close("CellTraxx Tuning Image"); exit("An unidentified error caused celltraxx.exe to stop running."); }

	tuning_image_code = 1; 											// First image is default
	button_number = 0; 
	while( button_number != 10 && button_number != 12 ){ 			// Running until [Run] or [Cancel] button pressed 

		setBatchMode(true); 										// To avoid flickering image on screen
		open(tun_dir+"celltraxx_tuning_image.bmp");
		run("Copy");
		close("celltraxx_tuning_image.bmp"); 	
		selectWindow("CellTraxx Tuning Image");
		run("Paste");		
		setBatchMode(false); 										// Needed to display tuning image on screen
		Roi.remove; 												// Removing the selection frame which would otherwise show in the image after paste
		getDimensions(img__width, img_height, channels, slices, frames);
		makeRectangle(left___crop_margin-1, inset_height-1+top____crop_margin, img__width-left___crop_margin-right__crop_margin+1, img_height-inset_height+1-top____crop_margin-bottom_crop_margin);
		if( selectionType > -1 ) {				// selectionType returns -1 if no selection exists 
			Roi.setStrokeColor("magenta");		// Setting the crop margin selection color and line width
			Roi.setStrokeWidth(1);
		}
		
		// Drawing numerical values right aligned into image with grey background 
		setJustification("right");
		gaussian_filter_radius_string = String.format("% 4.0f", gaussian_filter_radius);	setColor(  0,  0,  0);	drawString(gaussian_filter_radius_string, col_1_right, row_1_mid+font_vertical_shift, "#808080");	
		smallest_cell_diameter_string = String.format("% 3.0f", smallest_cell_diameter);	setColor(255,255,  0);	drawString(smallest_cell_diameter_string, col_1_right, row_2_mid+font_vertical_shift, "#808080");	
		largest__cell_diameter_string = String.format("% 4.0f", largest__cell_diameter);	setColor(255,  0,  0);	drawString(largest__cell_diameter_string, col_1_right, row_3_mid+font_vertical_shift, "#808080");	
		cutting__cell_diameter_string = String.format("% 4.0f", cutting__cell_diameter);	setColor(  0,  0,255);	drawString(cutting__cell_diameter_string, col_1_right, row_4_mid+font_vertical_shift, "#808080");	
	
		// Drawing a green frame around the curent tuning image button: First / Middle / Last
		setColor(  0,255,  0);	
		if(      tuning_image_code == 1 ) drawRect(col_3__left-1, tic_1_top-1, col_3_right-col_3__left+3, tic_1_bot-tic_1_top+3);	
		else if( tuning_image_code == 2 ) drawRect(col_3__left-1, tic_2_top-1, col_3_right-col_3__left+3, tic_2_bot-tic_2_top+3);	
		else if( tuning_image_code == 3 ) drawRect(col_3__left-1, tic_3_top-1, col_3_right-col_3__left+3, tic_3_bot-tic_3_top+3);	

/*		// Option to draw lines around button areas to check location on image 
		// + and - buttons
		drawRect(col_2__left, row_1_top, col_2_width, row_1_mid-row_1_top+1);	
		drawRect(col_2__left, row_1_mid, col_2_width, row_1_bot-row_1_mid+1);
		drawRect(col_2__left, row_2_top, col_2_width, row_2_mid-row_2_top+1);		
		drawRect(col_2__left, row_2_mid, col_2_width, row_2_bot-row_2_mid+1);
		drawRect(col_2__left, row_3_top, col_2_width, row_3_mid-row_3_top+1);	
		drawRect(col_2__left, row_3_mid, col_2_width, row_3_bot-row_3_mid+1);
		drawRect(col_2__left, row_4_top, col_2_width, row_4_mid-row_4_top+1);	
		drawRect(col_2__left, row_4_mid, col_2_width, row_4_bot-row_4_mid+1);
		// Update button
		drawRect(updat__left, row_1_top, updat_width, row_4_bot-row_1_top+1);
		// View image buttons
		drawRect(col_3__left, tic_1_top, col_3_right-col_3__left+1, tic_1_bot-tic_1_top+1);	
		drawRect(col_3__left, tic_2_top, col_3_right-col_3__left+1, tic_2_bot-tic_2_top+1);	
		drawRect(col_3__left, tic_3_top, col_3_right-col_3__left+1, tic_3_bot-tic_3_top+1);	
		// Run + Reset + Cancle buttons
		drawRect(col_4__left, row_1_top, col_4_right-col_4__left+1, row_1_bot-row_1_top+1);	
		drawRect(col_4__left, row_2_top, col_4_right-col_4__left+1, row_2_bot-row_2_top+1);	
		drawRect(col_4__left, row_3_top, col_4_right-col_4__left+1, row_3_bot-row_3_top+1);		
		drawRect(col_4__left, row_4_top, col_4_right-col_4__left+1, row_4_bot-row_4_top+1);		// */ 

		setLineWidth(1);
		button_number = 0;
		while( button_number == 0 ) {
			flags = 0;														// Making sure below while() loop is always entered at the start so the cursor and mouse are read initially
			while( flags&16==0 ) getCursorLoc(x, y, z, flags);				// Waiting for left click and then reading cursor coordinates (flags = 16 when Left Mouse Button is pressed)
			if( x>col_2__left && x<col_2_right && y>row_1_top && y<row_1_mid && gaussian_filter_radius<100                    ) { gaussian_filter_radius += 1; button_number =  4; }
			if( x>col_2__left && x<col_2_right && y>row_1_mid && y<row_1_bot && gaussian_filter_radius>2                      ) { gaussian_filter_radius -= 1; button_number =  5; }
			if( x>col_2__left && x<col_2_right && y>row_2_top && y<row_2_mid && smallest_cell_diameter<largest__cell_diameter ) { smallest_cell_diameter += 1; button_number =  6; }
			if( x>col_2__left && x<col_2_right && y>row_2_mid && y<row_2_bot && smallest_cell_diameter>2                      ) { smallest_cell_diameter -= 1; button_number =  7; }
			if( x>col_2__left && x<col_2_right && y>row_3_top && y<row_3_mid && largest__cell_diameter<200                    ) { largest__cell_diameter += 1; button_number =  8; }
			if( x>col_2__left && x<col_2_right && y>row_3_mid && y<row_3_bot && largest__cell_diameter>smallest_cell_diameter ) { largest__cell_diameter -= 1; button_number =  9; }
			if( x>col_2__left && x<col_2_right && y>row_4_top && y<row_4_mid && cutting__cell_diameter<200                    ) { cutting__cell_diameter += 1; button_number = 15; }
			if( x>col_2__left && x<col_2_right && y>row_4_mid && y<row_4_bot && cutting__cell_diameter>0                      ) { cutting__cell_diameter -= 1; button_number = 16; }
			if( x>col_3__left && x<col_3_right && y>tic_1_top && y<tic_1_bot ) button_number = tuning_image_code = 1;  	// First  image  
			if( x>col_3__left && x<col_3_right && y>tic_2_top && y<tic_2_bot ) button_number = tuning_image_code = 2;  	// Middle image
			if( x>col_3__left && x<col_3_right && y>tic_3_top && y<tic_3_bot ) button_number = tuning_image_code = 3;  	// Last   image
			if( x>col_4__left && x<col_4_right && y>row_1_top && y<row_1_bot ) button_number = 10;   // [Run]
			if( x>col_4__left && x<col_4_right && y>row_2_top && y<row_2_bot ) button_number = 11;   // [Reset]
			if( x>col_4__left && x<col_4_right && y>row_3_top && y<row_3_bot ) button_number = 12;   // [Cancel]
			if( x>col_4__left && x<col_4_right && y>row_4_top && y<row_4_bot ) button_number = 13;   // [Back]
			if( x>updat__left && x<updat_right && y>row_1_top && y<row_4_bot ) button_number = 14;   // [Update]

			// Reading selection coordinates in case the cropping frame is being moved in order to set new crop margins 
			if( selectionType > -1 ) {												// selectionType returns -1 if no selection exists 
				getSelectionBounds(roi_x, roi_y, roi__width, roi_height);
				if( roi__width > 33 && roi_height > 33 && roi_y > inset_height-1 ){			// Setting minimum accepted selection size somewhat arbitrarily
					left___crop_margin = roi_x + 1;	
					top____crop_margin = roi_y - inset_height + 1;
					right__crop_margin = img__width - roi__width - roi_x;
					bottom_crop_margin = img_height - roi_height - roi_y;
					// print("Image x y", img__width, img_height, "ROI x y", roi_x, roi_y, "ROI w h", roi__width, roi_height, "Crop margins l r t b", left___crop_margin, right__crop_margin, top____crop_margin, bottom_crop_margin);
				} 
			}     
		}
		// Drawing numerical values right aligned into image with grey background 
		gaussian_filter_radius_string = String.format("% 4.0f", gaussian_filter_radius);	setColor(  0,  0,  0);	drawString(gaussian_filter_radius_string, col_1_right, row_1_mid+font_vertical_shift, "#808080");	
		smallest_cell_diameter_string = String.format("% 3.0f", smallest_cell_diameter);	setColor(255,255,  0);	drawString(smallest_cell_diameter_string, col_1_right, row_2_mid+font_vertical_shift, "#808080");	
		largest__cell_diameter_string = String.format("% 4.0f", largest__cell_diameter);	setColor(255,  0,  0);	drawString(largest__cell_diameter_string, col_1_right, row_3_mid+font_vertical_shift, "#808080");	
		cutting__cell_diameter_string = String.format("% 4.0f", cutting__cell_diameter);	setColor(  0,  0,255);	drawString(cutting__cell_diameter_string, col_1_right, row_4_mid+font_vertical_shift, "#808080");	

		// Drawing a green frame around the curent tuning image button: First / Middle / Last, after first having drawn grey frames to remove the previous green frame 
		if( button_number > 0 && button_number < 4 ) { 
			setColor(128,128,128); 
			drawRect(col_3__left-1, tic_1_top-1, col_3_right-col_3__left+3, tic_1_bot-tic_1_top+3);	
			drawRect(col_3__left-1, tic_2_top-1, col_3_right-col_3__left+3, tic_2_bot-tic_2_top+3);	
			drawRect(col_3__left-1, tic_3_top-1, col_3_right-col_3__left+3, tic_3_bot-tic_3_top+3);	
			setColor(  0,255,  0);	
			if(      tuning_image_code == 1 ) drawRect(col_3__left-1, tic_1_top-1, col_3_right-col_3__left+3, tic_1_bot-tic_1_top+3);	
			else if( tuning_image_code == 2 ) drawRect(col_3__left-1, tic_2_top-1, col_3_right-col_3__left+3, tic_2_bot-tic_2_top+3);	
			else if( tuning_image_code == 3 ) drawRect(col_3__left-1, tic_3_top-1, col_3_right-col_3__left+3, tic_3_bot-tic_3_top+3);	
		}

		// Reversing the shaddowing of the button that was pressed to make it look like it moved inward
		if(      button_number ==  4 ) { setColor(64,64,64); drawLine(col_2__left, row_1_top  , col_2_right, row_1_top  ); drawLine(col_2__left, row_1_top  , col_2__left, row_1_mid-1); setColor(192,192,192); drawLine(col_2__left, row_1_mid-1, col_2_right, row_1_mid-1); drawLine(col_2_right, row_1_top+1, col_2_right, row_1_mid-1); }
		else if( button_number ==  5 ) { setColor(64,64,64); drawLine(col_2__left, row_1_mid+1, col_2_right, row_1_mid+1); drawLine(col_2__left, row_1_mid+1, col_2__left, row_1_bot  ); setColor(192,192,192); drawLine(col_2__left, row_1_bot  , col_2_right, row_1_bot  ); drawLine(col_2_right, row_1_mid+2, col_2_right, row_1_bot  ); }
		else if( button_number ==  6 ) { setColor(64,64,64); drawLine(col_2__left, row_2_top  , col_2_right, row_2_top  ); drawLine(col_2__left, row_2_top  , col_2__left, row_2_mid-1); setColor(192,192,192); drawLine(col_2__left, row_2_mid-1, col_2_right, row_2_mid-1); drawLine(col_2_right, row_2_top+1, col_2_right, row_2_mid-1); }
		else if( button_number ==  7 ) { setColor(64,64,64); drawLine(col_2__left, row_2_mid+1, col_2_right, row_2_mid+1); drawLine(col_2__left, row_2_mid+1, col_2__left, row_2_bot  ); setColor(192,192,192); drawLine(col_2__left, row_2_bot  , col_2_right, row_2_bot  ); drawLine(col_2_right, row_2_mid+2, col_2_right, row_2_bot  ); }
		else if( button_number ==  8 ) { setColor(64,64,64); drawLine(col_2__left, row_3_top  , col_2_right, row_3_top  ); drawLine(col_2__left, row_3_top  , col_2__left, row_3_mid-1); setColor(192,192,192); drawLine(col_2__left, row_3_mid-1, col_2_right, row_3_mid-1); drawLine(col_2_right, row_3_top+1, col_2_right, row_3_mid-1); }
		else if( button_number ==  9 ) { setColor(64,64,64); drawLine(col_2__left, row_3_mid+1, col_2_right, row_3_mid+1); drawLine(col_2__left, row_3_mid+1, col_2__left, row_3_bot  ); setColor(192,192,192); drawLine(col_2__left, row_3_bot  , col_2_right, row_3_bot  ); drawLine(col_2_right, row_3_mid+2, col_2_right, row_3_bot  ); }
		else if( button_number == 15 ) { setColor(64,64,64); drawLine(col_2__left, row_4_top  , col_2_right, row_4_top  ); drawLine(col_2__left, row_4_top  , col_2__left, row_4_mid-1); setColor(192,192,192); drawLine(col_2__left, row_4_mid-1, col_2_right, row_4_mid-1); drawLine(col_2_right, row_4_top+1, col_2_right, row_4_mid-1); }
		else if( button_number == 16 ) { setColor(64,64,64); drawLine(col_2__left, row_4_mid+1, col_2_right, row_4_mid+1); drawLine(col_2__left, row_4_mid+1, col_2__left, row_4_bot  ); setColor(192,192,192); drawLine(col_2__left, row_4_bot  , col_2_right, row_4_bot  ); drawLine(col_2_right, row_4_mid+2, col_2_right, row_4_bot  ); }
		else if( button_number == 10 ) { setColor(64,64,64); drawLine(col_4__left, row_1_top  , col_4_right, row_1_top  ); drawLine(col_4__left, row_1_top  , col_4__left, row_1_bot  ); setColor(192,192,192); drawLine(col_4__left, row_1_bot  , col_4_right, row_1_bot  ); drawLine(col_4_right, row_1_top+1, col_4_right, row_1_bot  ); wait(555); }
		else if( button_number == 11 ) { setColor(64,64,64); drawLine(col_4__left, row_2_top  , col_4_right, row_2_top  ); drawLine(col_4__left, row_2_top  , col_4__left, row_2_bot  ); setColor(192,192,192); drawLine(col_4__left, row_2_bot  , col_4_right, row_2_bot  ); drawLine(col_4_right, row_2_top+1, col_4_right, row_2_bot  ); wait(222); }
		else if( button_number == 12 ) { setColor(64,64,64); drawLine(col_4__left, row_3_top  , col_4_right, row_3_top  ); drawLine(col_4__left, row_3_top  , col_4__left, row_3_bot  ); setColor(192,192,192); drawLine(col_4__left, row_3_bot  , col_4_right, row_3_bot  ); drawLine(col_4_right, row_3_top+1, col_4_right, row_3_bot  ); wait(555); }
		else if( button_number == 13 ) { setColor(64,64,64); drawLine(col_4__left, row_4_top  , col_4_right, row_4_top  ); drawLine(col_4__left, row_4_top  , col_4__left, row_4_bot  ); setColor(192,192,192); drawLine(col_4__left, row_4_bot  , col_4_right, row_4_bot  ); drawLine(col_4_right, row_4_top+1, col_4_right, row_4_bot  ); wait(555); }
		else if( button_number == 14 ) { setColor(64,64,64); drawLine(updat__left, row_1_top  , updat_right, row_1_top  ); drawLine(updat__left, row_1_top  , updat__left, row_4_bot  ); setColor(192,192,192); drawLine(updat__left, row_4_bot  , updat_right, row_4_bot  ); drawLine(updat_right, row_1_top+1, updat_right, row_4_bot  ); wait(222); }
		close("\\Others"); 

		// Performing actions for the buttons which require updates to celltraxx_defaults.txt 
		if( button_number == 10 ) {			// [Run]				
			print("Finished interactive tuning. Running CellTraxx with the current settings.");
			perform_interactive_tuning = false; 				// Turning off tuning mode so CellTrax will do a full run this time 
			tuning_image_code = 0; 
		} else if( button_number == 11 ) { 	// [Reset]
			print("The tuning process was reset. Using values from last run."); 
			gaussian_filter_radius = original_gaussian_filter_radius;	
			smallest_cell_diameter = original_smallest_cell_diameter;
			largest__cell_diameter = original_largest__cell_diameter;
			cutting__cell_diameter = original_cutting__cell_diameter;
			top____crop_margin     = original_top____crop_margin;
			bottom_crop_margin     = original_bottom_crop_margin;
			left___crop_margin     = original_left___crop_margin;
			right__crop_margin     = original_right__crop_margin;
			// Drawing numerical values right aligned into image with grey background 
			gaussian_filter_radius_string = String.format("% 4.0f", gaussian_filter_radius);	setColor(  0,  0,  0);	drawString(gaussian_filter_radius_string, col_1_right, row_1_mid+font_vertical_shift, "#808080");	
			smallest_cell_diameter_string = String.format("% 3.0f", smallest_cell_diameter);	setColor(255,255,  0);	drawString(smallest_cell_diameter_string, col_1_right, row_2_mid+font_vertical_shift, "#808080");	
			largest__cell_diameter_string = String.format("% 4.0f", largest__cell_diameter);	setColor(255,  0,  0);	drawString(largest__cell_diameter_string, col_1_right, row_3_mid+font_vertical_shift, "#808080");	
			cutting__cell_diameter_string = String.format("% 4.0f", cutting__cell_diameter);	setColor(  0,  0,255);	drawString(cutting__cell_diameter_string, col_1_right, row_4_mid+font_vertical_shift, "#808080");	
		} else if( button_number == 13 ) { 	// [Back] to input dialog
			close("*"); 
			Write_Settings_To_File(); 							// Saving currents settings before restarting 
			runMacro(sys_dir+"celltraxx_imagej_macro.ijm"); 	// For simplicity, just running macro recursively when [Back] button is pressed 
			exit; 
		} 

		// Updating name of results folder based on current settings and generating full name
		default_part_of_results_folder_name = "GFR="+gaussian_filter_radius+"um_d="+smallest_cell_diameter+"-"+largest__cell_diameter+"um_cut="+cutting__cell_diameter+"um"; 
		if( first_part_of_results_folder_name == "" ) { results_folder_name = default_part_of_results_folder_name+ffc_string;    first_part_of_results_folder_name = "#"; }
		else 							   			    results_folder_name = first_part_of_results_folder_name+"_"+default_part_of_results_folder_name+ffc_string;

		Write_Settings_To_File(); 		// Writing settings back to file to store any changes 

		if( button_number == 12 ) { 	// [Cancel]
			close("*"); 
			exit;  
		} else if( button_number == 14 || (button_number > 0 && button_number < 4 ) ) {						// [Update], [First], [Middle] or [Last] pressed
			exec("cmd /c start "+sys_dir+"celltraxx");    													// "start" will open DOS window for output during execution 
			filestring = File.openAsString(sys_dir+"celltraxx_error_messages.txt");							// Checking if celltraxx.exe left a warning or error message, if so exiting and writing the error text to a message box
			if(      lengthOf(filestring) > 11 ) { close("CellTraxx Tuning Image"); exit(filestring); }  	// Exiting macro after closing tuning window and writing the error message from celltraxx.exe to message box
			else if( lengthOf(filestring) == 0 ) { close("CellTraxx Tuning Image"); exit("An unidentified error caused celltraxx.exe to stop running."); }
		}
		while( flags&16 > 0 ) getCursorLoc(x, y, z, flags);													// Waiting until left mouse button is lifted before continuing
	}
	close(); 
}

// CALLING THE CELLTRAXX EXCECUTABLE AND LETTING IT PROCESS ALL VIDEOS UNTIL FINISHED 
print("\nRunning the CellTraxx program. This can take some time. Check status in the command window C:\\celltraxx_system\\celltraxx.exe");
exec("cmd /c start "+sys_dir+"/celltraxx");   														// "start" will open command window for output during execution 
filestring = File.openAsString(sys_dir+"celltraxx_error_messages.txt");								// Checking if celltraxx.exe left a warning or error message, if so exiting and writing the error text to a message box
	if(      lengthOf(filestring) > 11 ) { close("CellTraxx Tuning Image"); exit(filestring); }  	// Exiting macro after closing tuning window and writing the error message from celltraxx.exe to message box
	else if( lengthOf(filestring) == 0 ) { close("CellTraxx Tuning Image"); exit("An unidentified error caused celltraxx.exe to stop running."); }

setBatchMode(true);

// POST PROCESSING DATA FROM CELLTRAXX.EXE: IMAGE SERIES MADE INTO VIDEOS AND SELECTED .CSV RESULTS FILES MADE INTO CHARTS AND SAVED AS .PNG IMAGES
new_dir = dat_dir + results_folder_name + "/";
print("\nPlacing results in folder " + new_dir); 
filestring = File.openAsString(new_dir+"celltraxx_bitmaps2process.txt");
textlines = split(filestring, "\n");

print("\nGenerating .avi files:"); 
i=1;
while( textlines[i]!="stop" ){
	
	columns=split(textlines[i],"\t");		// Splitting text line into name, first image and last image 
	rootname    =          columns[0];
	first_image = parseInt(columns[1]);
	last__image = parseInt(columns[2]);
 	print("First Image :", first_image, " Increment :", increment, " Last Image :", last__image); 
	
 	if( make_identified_cell_videos ) {
		for (j = first_image; j <= last__image; j += increment){
			numberstring = IJ.pad(j, 4);
			bmpname = new_dir + rootname + "_05_" + numberstring + "_Identified.bmp"; 
			open(bmpname); 
			if( !write_identified_cell_images ) dummy = File.delete(bmpname);
		}
		run("Images to Stack", "name=[] title=[]");
		aviname = new_dir + rootname + "_Video_Identified.avi"; 	
		save_code = "compression=None frame=4 save=" + aviname; 	
		run("AVI... ", save_code);
		print("Saved AVI video with identified cells as ", aviname); 
		close(); 
 	}
	if( make_matched_cell_videos ) {
		for (j = first_image; j <= last__image; j += increment){
			numberstring = IJ.pad(j, 4);
			bmpname = new_dir + rootname + "_06_" + numberstring + "_Matched.bmp"; 
			open(bmpname);
			if( !write_matched_cell_images ) dummy = File.delete(bmpname);
		}
		run("Images to Stack", "name=[] title=[]");
		aviname = new_dir + rootname + "_Video_Matched.avi"; 	
		save_code = "compression=None frame=4 save=" + aviname; 	
		run("AVI... ", save_code);
		print("Saved AVI video with matched cells as ", aviname); 
		close(); 
 	}
 	if( make_valid_track_videos ) {			
		for (j = first_image; j <= last__image; j += increment){
			numberstring = IJ.pad(j, 4);
			bmpname = new_dir + rootname + "_07_" + numberstring + "_Valid_Track.bmp"; 
			open(bmpname);
			if( !write_valid_track_images ) dummy = File.delete(bmpname);
		}
		run("Images to Stack", "name=[] title=[]");
		aviname = new_dir + rootname + "_Video_Valid_Track.avi"; 	
		save_code = "compression=None frame=4 save=" + aviname; 	
		run("AVI... ", save_code);
		print("Saved AVI video with valid track cells as ", aviname); 
		close(); 
 	}

	i++;

	// ======== VELOCITY HISTOGRAM ============
	csvstring = File.openAsString(new_dir + rootname + "_Velocity-histogram.csv");
	csvlines = split(csvstring, "\n");
	x_values = newArray(nof_bins_in_histogram);
	y_values = newArray(nof_bins_in_histogram);
	p_values = newArray(nof_bins_in_histogram);
	y_max = 0;
	for( n = 1; n<nof_bins_in_histogram; n++ ){
		entries=split(csvlines[n],",");
		x_values[n] = parseFloat(entries[0]);
		y_values[n] = parseFloat(entries[1]);
		p_values[n] = parseFloat(entries[2]);
		if( y_values[n] > y_max ) y_max = y_values[n];
	}	
	y_max = 1000*( floor(y_max/1000)+1); // Rounding upwards to nearest 1000 counts
	Plot.create("Velocity Histogram", "Cell velocity [um/min]", "Counts");
	Plot.setLineWidth(2);
	Plot.setFontSize(30);	// setFontSize(40, "bold");
	Plot.setLimits(0, highest_cell_velocity, 0, y_max);
	Plot.setFrameSize(800,600);
	Plot.setColor("red" );		Plot.add("bars", x_values, y_values);
	Plot.setColor("blue");		Plot.add("line", x_values, p_values);
	Plot.show();
	save(new_dir + rootname + "_Velocity-histogram.png"); 
	close();

	// ======== VELOCITY VERSUS TIME ============
	csvstring = File.openAsString(new_dir + rootname + "_Velocity-matrix.csv");
	csvlines = split(csvstring, "\n");
	x_values = newArray(csvlines.length-9);		// There are 9 lines of header
	y_values = newArray(csvlines.length-9);
	x_max = y_max = 0;
	n = 9;		// Skipping  header lines of Velocity-matrix.csv file 
	m = 0; 
	while( n < csvlines.length ){	
		entries=split(csvlines[n],",");
		x_values[m] = parseFloat(entries[2]);
		y_values[m] = parseFloat(entries[5]);
		if( n == 9 )              x_min = x_values[m];
		if( x_values[m] > x_max ) x_max = x_values[m];
		if( y_values[m] > y_max ) y_max = y_values[m];
		n++;
		m++;
	}	
	x_min = (floor(x_min)  ); 	// Rounding   down  to nearest whole hour
	x_max = (floor(x_max)+1); 	// Rounding upwards to nearest whole hour
	y_max = (floor(y_max)+1); 	// Rounding upwards to nearest whole µm/min
//	x_min = 0; x_max = 10; y_max = 2; 			//  TO OVERRIDE THE AUTOMATIC AXIS LIMITS, ADD CONSTANT VALUES HERE AND SIMILARLY FOR THE OTHER PLOTS 
	Plot.create("Velocity versus Time", "Time [h]", "Cell velocity [um/min]");
	Plot.setLineWidth(2);
	Plot.setFontSize(30); 
	Plot.setLimits(x_min, x_max, 0, y_max);
	Plot.setFrameSize(800,600);
	Plot.setColor("red" );
	Plot.add("line", x_values, y_values);
	Plot.show();
	save(new_dir + rootname + "_Velocity-vs-Time.png"); 
	close();
	
	// ======== POSITION DATA plotted from COMMON ORIGIN ============
	csvstring = File.openAsString(new_dir + rootname + "_Position-data.csv");
	csvlines = split(csvstring, "\n");
	max_nof_images_in_video = floor( (last__image-first_image+1)/increment + 2 );
	x_pos = newArray(max_nof_images_in_video);
	y_pos = newArray(max_nof_images_in_video);
	abs_max = 0;									// Going through all the data once to find the largest absolute value in x and y to allow setting correct plot size 
	last_x_in_track = newArray(3000);				// Dimension should be equal to MAXCELLS in celltraxx.c
	last_y_in_track = newArray(3000);
	hex_color_track = newArray(3000);
	m = 1; 
	while( m < csvlines.length ){
		if( csvlines[m] != "" ) { 
			entries = split(csvlines[m],",");
			x = abs( parseFloat(entries[7]) );
			y = abs( parseFloat(entries[8]) );
			if( x > abs_max ) abs_max = x;
			if( y > abs_max ) abs_max = y;
		} 
		m++;
	}
	abs_max = 100*( floor(abs_max/100)+1);  		// Rounding upwards to nearest 100 µm
	Plot.create("Cell tracks plotted from common origin", "Horizontal position X [um]", "Vertical position Y [um]\n");
	Plot.setLineWidth(2);
	Plot.setFontSize(30);
	Plot.setLimits(-abs_max, abs_max, -abs_max, abs_max);
	Plot.setFrameSize(800,800);
	m = 1; 											// Going through all the data again and storing x and y values in arrays for plotting of all the tracks separately
	n = 0;
	while( m < csvlines.length ){
		if( csvlines[m] != "" ) { 
			entries = split(csvlines[m],",");
			x_pos[n] = parseFloat(entries[7]);
			y_pos[n] = parseFloat(entries[8]);
			n++;
		} else if( n > 0 ) {									// Needed this test in cases of two blank lines at the end of the file to avoid index [-1] below
			for( k = n; k < max_nof_images_in_video; k++ ) {	// Filling rest of array with the trace line end coordinates to avoid lines returning to the origin in the plots
				x_pos[k] = x_pos[n-1];
				y_pos[k] = y_pos[n-1];
			}
			Plot.setColor("red");
			Plot.add("line", x_pos, y_pos);
			n = 0;
		}
		m++;
	}
	Plot.show();
	save(new_dir + rootname + "_Position-data_Common-origin.png"); 
	close();
	// Making a new plot where the track lines have the same colors as in the .avi videos and then adding track numbers for the long tracks 
	Plot.create("Cell tracks plotted from common origin", "Horizontal position X [um]", "Vertical position Y [um]\n");
	Plot.setLineWidth(2);
	Plot.setFontSize(30);
	Plot.setLimits(-abs_max, abs_max, -abs_max, abs_max);
	Plot.setFrameSize(800,800);
	m = 1; 											// Going through all the data again and storing x and y values in arrays for plotting of all the tracks separately
	n = 0;
	nof_tracks = 0;
	while( m < csvlines.length ){
		if( csvlines[m] != "" ) { 
			entries = split(csvlines[m],",");
			x_pos[n] = parseFloat(entries[7]);
			y_pos[n] = parseFloat(entries[8]);
			hex_color = '#'+entries[14];						// To allow plotting the lines with the same colors as in videos
			n++;
		} else if( n > 0 ) {									// Needed this test in cases of two blank lines at the end of the file to avoid index [-1] below
			for( k = n; k < max_nof_images_in_video; k++ ) {	// Filling rest of array with the trace line end coordinates to avoid lines returning to the origin in the plots
				x_pos[k] = x_pos[n-1];
				y_pos[k] = y_pos[n-1];
			}
			Plot.setColor(hex_color);
			Plot.add("line", x_pos, y_pos);
			last_x_in_track[nof_tracks]=x_pos[n-1];				// Storing track end point coordinates in separate arrays for adding track numbers below
			last_y_in_track[nof_tracks]=y_pos[n-1];
			hex_color_track[nof_tracks]=hex_color;
			nof_tracks++;		
			n = 0;
		}
		m++;
	}
	// Adding track numbers with correct colors at the end of long tracks. Using trick with the Plot.add("code...") feature to get the right colors
	x_label = newArray(1);
	y_label = newArray(1);
	for( n=0; n<nof_tracks; n++ ) {
		if( sqrt( last_x_in_track[n]*last_x_in_track[n] + last_y_in_track[n]*last_y_in_track[n] ) > 0.4*abs_max ) {		// Adding number only if track end point is > 40 % of plot frame half-size away from the origin
			x_label[0] = last_x_in_track[n];
			y_label[0] = last_y_in_track[n];
		    Plot.setColor(hex_color_track[n]);
		    Plot.add("code: setFont('sanserif',20,'bold');drawString(''+"+d2s(n+1,0)+",x-8,y+8);", x_label, y_label); 
		}
	}
	Plot.show();
	save(new_dir + rootname + "_Position-data_Common-origin_Long-tracks-numbered.png"); 
	close();
	
	// ======== POSITION DATA plotted as IN IMAGE ============
	csvstring = File.openAsString(new_dir + rootname + "_Position-data.csv");
	csvlines = split(csvstring, "\n");
	m = 1; 
	while( m < csvlines.length ){
		if( csvlines[m] != "" ) { 
			entries = split(csvlines[m],",");
			x = parseFloat(entries[5]);
			y = parseFloat(entries[6]);
			if( x > x_max ) x_max = x;
			if( y > y_max ) y_max = y;
		} 
		m++;
	}
	x_max = 100*( floor(x_max/100)+1);  		// Rounding upwards to nearest 100 µm
	y_max = 100*( floor(y_max/100)+1);  
	Plot.create("Cell tracks plotted as in image and numbered as in Position-data.csv file", "Horizontal position X [um]", "Vertical position Y [um]\n");
	Plot.setLineWidth(2);
	Plot.setFontSize(30);
	Plot.setLimits(0, x_max, 0, y_max);
	Plot.setFrameSize(x_max,y_max);
	Plot.setBackgroundColor("#999999")
	m = 1; 											// Going through all the data again and storing x and y values in arrays for plotting of all the tracks separately
	n = 0;
	nof_tracks = 0;
	while( m < csvlines.length ){
		if( csvlines[m] != "" ) { 
			entries = split(csvlines[m],",");
			x_pos[n] = parseFloat(entries[5]);
			y_pos[n] = parseFloat(entries[6]);
//			hex_color = '#'+entries[14];						// To allow plotting the lines with the same colors as in videos
			n++;
		} else if( n > 0 ) {									// Needed this test in cases of two blank lines at the end of the file to avoid index [-1] below
			for( k = n; k < max_nof_images_in_video; k++ ) {	// Filling rest of array with the trace line end coordinates to avoid lines returning to the origin in the plots
				x_pos[k] = x_pos[n-1];
				y_pos[k] = y_pos[n-1];
			}
			Plot.setColor(hex_color_track[nof_tracks]);
			Plot.add("line", x_pos, y_pos);
			last_x_in_track[nof_tracks]=x_pos[n-1];				// Storing track end point coordinates in separate arrays for adding track numbers below
			last_y_in_track[nof_tracks]=y_pos[n-1];
			nof_tracks++;		
			n = 0;
		}
		m++;
	}
	Plot.show();	
	save(new_dir + rootname + "_Position-data-as-in-image.png"); 
	// Adding track numbers with correct colors at the end of each track. Using trick with the Plot.add("code...") feature to get the right colors
	for( n=0; n<nof_tracks; n++ ) {
		x_label[0] = last_x_in_track[n];
		y_label[0] = last_y_in_track[n];
	    Plot.setColor(hex_color_track[n]);
	    Plot.add("code: setFont('sanserif',20,'bold');drawString(''+"+d2s(n+1,0)+",x-8,y+8);", x_label, y_label); 
	}
	save(new_dir + rootname + "_Position-data-as-in-image_Tracks-numbered.png"); 
	close();

	// ======== DIAMETER HISTOGRAM ============
	csvstring = File.openAsString(new_dir + rootname + "_Diameter-histogram.csv");
	csvlines = split(csvstring, "\n");
	d_values = newArray(nof_bins_in_histogram);
	c_values = newArray(nof_bins_in_histogram);
	a_values = newArray(nof_bins_in_histogram);
	y_max = 0;
	for( n = 1; n<nof_bins_in_histogram; n++ ){
		entries=split(csvlines[n],",");
		d_values[n] = parseFloat(entries[0]);	// Diameter
		c_values[n] = parseFloat(entries[1]);	// Candidate cells
		a_values[n] = parseFloat(entries[2]);	// Identified  cells
		if( c_values[n] > y_max ) y_max = c_values[n];
		if( a_values[n] > y_max ) y_max = a_values[n];
	}	
	y_max = 200*( floor(y_max/200)+1); // Rounding upwards to nearest 200 counts
	Plot.create("Diameter Histogram", "Equivalent cell diameter [um]", "Counts");
	Plot.setLineWidth(2);
	Plot.setFontSize(30);
	Plot.setLimits(0, largest__cell_diameter, 0, y_max);
	Plot.setFrameSize(800,600);
	Plot.setLineWidth(2);	Plot.setColor("#999999" );	Plot.add("line", d_values, c_values);
	Plot.setLineWidth(2);	Plot.setColor("#00BB00");	Plot.add("bars", d_values, a_values);
	Plot.setLegend("Candidate cells\tIdentified cells", "top-right");
	Plot.show();
	save(new_dir + rootname + "_Diameter-histogram.png"); 
	close();
		
}
setBatchMode(false);
print("\nFinished processing all .avi videos in the folder "+dat_dir+"\n");


// Function which writes the current settings to the file celltraxx_defaults.txt in the systems directory

function Write_Settings_To_File()
{
	//	print("Saving data to file celltraxx_defaults.txt"); 
	if( wound_healing_mode              ) string = "yes"; else string = "no"; File.saveString("Wound healing mode              "+string+"\n", sys_dir+"celltraxx_defaults.txt");		
	if( perform_flat_field_correction   ) string = "yes"; else string = "no"; File.append(    "Perform flat field correction   "+string,      sys_dir+"celltraxx_defaults.txt");		
	if( perform_image_shift_correction  ) string = "yes"; else string = "no"; File.append(    "Perform image shift correction  "+string,      sys_dir+"celltraxx_defaults.txt");		
	if( perform_interactive_tuning 	    ) string = "yes"; else string = "no"; File.append(    "Perform interactive tuning      "+string,      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Version                         "+version, 			          												 	      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Track smoothing iterations      "+track_smoothing_iterations,       												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Results folder drive letter     "+dat_drive_letter,                  												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Results folder name             "+results_folder_name,               												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "First part of folder name       "+first_part_of_results_folder_name, 												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Pixel size [um]                 "+pixel_size, 					      												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Gaussian filter radius [um]     "+gaussian_filter_radius, 		      												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Smallest cell diameter [um]     "+smallest_cell_diameter, 		      												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Largest  cell diameter [um]     "+largest__cell_diameter, 		      												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Cutting cell diameter [um]      "+cutting__cell_diameter, 		      												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Top    crop margin [pixels]     "+top____crop_margin, 			      												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Bottom crop margin [pixels]     "+bottom_crop_margin, 			      												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Left   crop margin [pixels]     "+left___crop_margin, 		      	  												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Right  crop margin [pixels]     "+right__crop_margin, 			      												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Time between images      [min]  "+time_between_images, 		      												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Highest cell velocity [um/min]  "+highest_cell_velocity, 		      												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Shortest cell track [images]    "+shortest_cell_track, 		      												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "First image number              "+first_image,                       												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Last image number               "+last_image,                        												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Image number increment          "+increment,                         												      sys_dir+"celltraxx_defaults.txt");		
	if( make_identified_cell_videos 	) string = "yes"; else string = "no"; File.append(    "Make identified cell videos     "+string,      sys_dir+"celltraxx_defaults.txt");		
	if( make_matched_cell_videos 		) string = "yes"; else string = "no"; File.append(    "Make matched cell videos        "+string,      sys_dir+"celltraxx_defaults.txt");		
	if( make_valid_track_videos 		) string = "yes"; else string = "no"; File.append(    "Make valid track videos         "+string,      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Tracking dot diameter [um]      "+tracking__dot_diameter, 		      												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Valid track images contrast     "+valid_track_image_contrast,        												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Scalebar color                  "+scalebar_color,                    												      sys_dir+"celltraxx_defaults.txt");		
	if( draw_cell_outline 	            ) string = "yes"; else string = "no"; File.append(    "Draw cell outline               "+string,      sys_dir+"celltraxx_defaults.txt");		
	if( draw_cell_track_line 			) string = "yes"; else string = "no"; File.append(    "Draw cell track line            "+string,      sys_dir+"celltraxx_defaults.txt");		
	if( write_mirror_margin_images 		) string = "yes"; else string = "no"; File.append(    "Write mirror margin images      "+string,      sys_dir+"celltraxx_defaults.txt");		
	if( write_shifted_images 			) string = "yes"; else string = "no"; File.append(    "Write shifted images            "+string,      sys_dir+"celltraxx_defaults.txt");		
	if( write_gaussian_smoothed_images 	) string = "yes"; else string = "no"; File.append(    "Write gaussian smoothed images  "+string,      sys_dir+"celltraxx_defaults.txt");		
	if( write_segmented_cell_images 	) string = "yes"; else string = "no"; File.append(    "Write segmented cell images     "+string,      sys_dir+"celltraxx_defaults.txt");		
	if( write_cut_cell_images 		    ) string = "yes"; else string = "no"; File.append(    "Write cut cell images           "+string,      sys_dir+"celltraxx_defaults.txt");		
	if( write_identified_cell_images 	) string = "yes"; else string = "no"; File.append(    "Write identified cell images    "+string,      sys_dir+"celltraxx_defaults.txt");		
	if( write_matched_cell_images 		) string = "yes"; else string = "no"; File.append(    "Write matched cell images       "+string,      sys_dir+"celltraxx_defaults.txt");		
	if( write_valid_track_images 		) string = "yes"; else string = "no"; File.append(    "Write valid track images        "+string,      sys_dir+"celltraxx_defaults.txt");		
	if( keep_cells_from_previous_image  ) string = "yes"; else string = "no"; File.append(    "Keep cells from previous image  "+string,      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Segmentation limit [SDs]        "+segmentation_limit, 			      												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "# bins in histogram             "+nof_bins_in_histogram, 		      												      sys_dir+"celltraxx_defaults.txt");		
	File.append(    "Tuning image code               "+tuning_image_code,                 												      sys_dir+"celltraxx_defaults.txt");		
}


