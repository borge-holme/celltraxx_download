#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "direct.h"					// Needed for creating directory using _mkdir()

#define	VERSION 4.6					// New in 4.6: Adjusted spaces in celltraxx_defaults.txt so it looks nicer. // New in 4.5: Drawing cutting lines in tuning image since macro has option to tune the value. // New in 4.4: Changed name for image numbers in command file, copying previous segmented image after scrapping large and small cells. // New in 4.3: New order of celltraxx_defaults.txt to match input window in ImageJ. No longer need to specify max length of image shift. // New in 4.2: All code merged into one where only the used routines from bmpio.c are included. // New in 4.0: Reading common imange number increment from celltraxx_defaults.txt and omitting it from celltraxx_bitmaps2process.txt // New in 3.9: Improved cutting function, temporary tuning folder. // New in 3.8: Funcion cut_concave_region() which tries to divide large regions by cutting at narrow convex points. New in 3.6: 1) New format for _Overall_summary.csv. 2) The Identified image is drawn with just the cell outline on top of the original greyscale image. 3) Text in celltraxx_defaults.txt changed to "Draw cell outline". // Version 3.5 = 3.4 but new number to match new layout in ImageJ tuning image. // New in 3.4 is flat field correction done in ImageJ macro but written in command file. // New in 3.3 is two separate folders C:\celltraxx_system and X:\celltraxx_data. Also writing error messages to file so ImageJ macro can read and display. // New in 3.2 is to search at all eight positions ±3 pixels away in image shift and connecting tracks which lack cell in only one image. New in 3.1 is writing both full name of results folder and also added name at bottom of celltraxx_default.txt. New in 3.0: Generating tuning image and exiting if tuning mode set by user. New in 2.1: Output of size histogram to allow easier setting of smallest and largest cell diameter. // New in 2.0: Generating the ..._DOS_Window_Table.csv for easy plotting of grey level and cell number data in Excel. Corrected scrapping of large regions (and eliminating erosion and dilation). User defined color dot diameter. Automatic and robust grey level segmentation limit by taking darkest level in background for first image in video. // New in 1.9 is option to smoothen cell positions through time in order to reduce back-and-forth movements of cells. This gives lower and more correct velocity values for stationary cells. // New in 1.8 is fractional image shift and true interpolation when shifting. // New in 1.7 is that the Segmentation_Limit is given directly in the command file (and not multiplied by the SD as in older versions). // New in 1.5 is that bottom of command file gives name of directory in which to place images and .csv results files. Default directory is "Results". // New in 1.4 is option to send results to separate folder and to run only a subset of original bitmaps with chosen image step. // New in 1.2 is top crop margin which is useful if images have black lines at both top and bottom. // New in 1.1 is reading initial bitmap filename and number of images from text file celltraxx_bitmaps2process.txt
#define MAXX      2400				// Max number of horizontal pixels in images (+ a little margin for filtering ). Using same size as defined for bitmaps in celltraxx_bmpio.c. 
#define MAXY      2400				// Max number of  vertical  pixels in images (+ a little margin for filtering )
#define MAXBINS    501				// Max number of bins in histogram
#define MAXLEVELS  500				// Max number of images in any series (stack) to process	
#define MAXCELLS  6000				// Maximum number of cells anticipated in one image. With 3000 cells, celltraxx.exe uses about 500 MB RAM, with 6000 cells it uses 850 MB. 				
#define MATCHED_CELLS 6*MAXCELLS	// Max number of individual cells the program can track through all the image levels 
#define SHIFTXY	100					// Number of pixels from origin to shift images when placing them into Expanded_Image[][], to allow easy filtering and sideways adjustment (if images not well aligned) 
#define	HUESPERCOLOR 20				// Number of different hues per main color (blue, cyan, green, yellow, red, magenta) to generate for marking of cells 
#define DIVISION_MARGIN  2			// The fraction of mismatch allowed when comparing the distances and areas of two new cells suspected of being the result of a cell division 
#define DIVISION_ANGLE 120			// The maximum allowed divergence angle between a vector from one of the new cells to its twin and to its mother. This ensures that the three cells lie nearly on a straight line.
#define TOO_SMALL_CELL   1			// Used for coding cells in tuning image 
#define TOO_LARGE_CELL   2			// Used for coding cells in tuning image 
#define TOO_FAR_FROM_WH_EDGE 3		// Used for coding cells in tuning image that are too far from the frontmost cell during wound healing 
#define CIRCLE_FIT_POINT_DIST 8		// Used in cut_concave_region() as nuber of points to each side of a concave minimum to use when fitting a circle to get local radius of curvature
#define CUT_MARGIN 3				// Distance of overlap in Edge_Point_Distances[][] to allow simple search for minima without edge effects.
#define CUT_MAX_EDGE_POINTS 3*MAXX	// Size of Edge_Point_X_px[] and Edge_Point_Y_px[] to hold coordinates of cutting edge points. Allowing a convoluted region with border longer than the side edges. Used in cut_concave_region().
#define MAX_SEARCH_SPIRAL_EDGE 9	// The longest edge of the "square shaped search spiral" to try before giving up the search for the next edge point in cut_concave_region()
#define WARNING 0					// Used in write_error_message() to distinguish warning and error
#define ERROR	1
#define TRUE    1
#define FALSE   0
#define STRING_LENGTH 200			// Max number of characters in text strings like file names, etc. 
#define CUT  -32762					// Special code for cut lines made by cut_concave_region() to allow display in images for convenience 
#define CUT2 -32763
#define FIT  -32764
#define FIT2 -32765
#define FIT3 -32766
#define BIG   32766					// The largest value for a short - 1
#define FILL_TUN(x,y) if( x>=0 && x<Bimage->xsize && y>=0 && y<Bimage->ysize &&     Tuning_Image[SHIFTXY+x][SHIFTXY+y]==fillnumber) {     Tuning_Image[SHIFTXY+x][SHIFTXY+y]=withnumber; count++; PointList_x[plistwrite]=x;	PointList_y[plistwrite]=y; plistwrite = (plistwrite+1)%(2*MAXY); }		// Used in fill_tuning()  
#define FILL_SEG(x,y) if( x>=0 && x<Bimage->xsize && y>=0 && y<Bimage->ysize &&  Segmented_Image[SHIFTXY+x][SHIFTXY+y]==fillnumber) {  Segmented_Image[SHIFTXY+x][SHIFTXY+y]=withnumber; count++; PointList_x[plistwrite]=x;	PointList_y[plistwrite]=y; plistwrite = (plistwrite+1)%(2*MAXY); }		// Used in fill_segmented()  

//*************************************** Definition of structures **********************************

typedef struct cell {
	int cellnumber;							// Number which identifies each physical cell through the various levels. Given to cells identified in first level and then to cells matched in subsequent levels. 							
	int area_px;							// Number of pixels in segmented area								
	int cx, cy;								// Center coordinates (in hole pixels)									
	int cm_x_first, cm_y_first;				// The first pixels in the cell which was encountered when searching through segmented cells. Needed when changing color during matching since center of mass is not always inside the colored area 								
	int min_x, max_x, min_y, max_y;			// The bounding box coordinates returned by the fill_segmented_box() function. Used for estimating circularity. 
	double dx, dy;							// Center coordinates (in pixels) but with decimals, used for calculating exact realphi and realtheta from the shift in position
	double fx, fy;							// Center coordinates (in pixels) used for temporary storage during time filtering of positions to get more correct velocities for slow moving cell
	double aspect_ratio;					// The height of the bounding box divided by the width. Will be close to 1 for a circular cell just before and after division. 
	double area_fraction;					// The area (area_px) of the cell divided by Pi/4 the area of the bounding box. Will be close to 1 if circular (or elliptical) cell and less if more elongated diagonally or irregular. 
} CELL;

typedef struct bmap {						// Defines a bitmap image containing:			
	char  path[STRING_LENGTH];				// Bitmap path									
	char  filename[STRING_LENGTH];			// Bitmap filename -extension					
	char  ext[9];							// Bitmap extension bimp						
	char  b, m;								// Should be "BM"								
	long int filesize;						// Number of bytes in file						
	long int reserved;						// =0											
	long int dataoffset;					// File offset to raster data					
	long int infoheadersize;				// =40											
	long int xsize, ysize;					// Number of pixels in x and y direction		
	long int xsize4, ysize4;				// Number of pixels in x and y direction		
	short int planes;						// Number of planes = 1							
	short int bitcount;						// =8 for 8 bit palletized 256 colors			
	long  int compression;					// =0 for no compression						
	long  int imagesize;					// Size of compressed image						
	long  int xres, yres;					// # pixels/m in x and y directions				
	long  int colorsused;					// # colors used								
	long  int colorsimport;					// # important colors							
	long  int orientation;					// From tiff images. First pixel is: 1=top left, 2= top right, 3=bottom right, 4=bottom left 
	long  int samplesperpixel;				// From tiff images. Number of channels in each pixel 
	unsigned char colortable[256][3];		// Translating color number to Blue Green Red	
	unsigned char image[MAXX][MAXY][3];		// Raster of image data	with last index Blue Green Red where image[0][0][] is lower left corner				
} BITMAP;

//**********************************  Global varialbes  *******************************************

char  BMPfilename[STRING_LENGTH];							// Full name of first image to read
char  ResultsFolderName[STRING_LENGTH];						// Name of (sub) folder in which to write all results from the run using the current command file 
char  InputFolderName[STRING_LENGTH];						// Name of type X:/celltraxx where the drive letter X is read from the celltraxx_defaults.txt to allow image data read from external disks if no room on local drive 
char  PreviousFirstImageName[STRING_LENGTH*2];				// Name of first image to be processed in first video, stored in order to check if it is necessary to find a new Segmentation_Grey_Limit_From_First_Image
char  Draw_CM, Draw_Outline, Draw_Track;					// Flags set in command file to determine what is included in the output images. 
char  Keep_Previous_Cells;									// Flags set in command file if using identified cells in the previous image to mark positions of candidate cells in current image. 
char  Perform_Tuning, Perform_FFC, Perform_Shifting;		// Flags set in command file for interactive tuning on first video, flat field correction and image shift correction. Since flat field correction is conveniently done in the CellTraxx.ijm macro before the avi video is converted to bitmaps, Perform_FFC does not have a function here. 
char  Use_SD_Based_Segmentation_Grey_Limit;					// Flag for overruling the auto segmentation limit function and instead using the grey mean value minus the user given Segmentation Limit value * the standard deviation (SD) of the image. Flag turned on by giving a negative Segmentation Limit 
char  Scalebar_Color[40], Add_Scalebar, Add_100um_image;	// Text from command file with desired color of scale bar, flags used for adding scale bar and text "100 µm" 
char  Write_0_MM_img, Write_1_SH_img;  
char  Write_2_SM_img, Write_3_SE_img, Write_4_CU_img;		// Flags set in command file to control output images written to file from various stages of the image processing. Normally enough with all except the last one. 
char  Write_5_ID_img, Write_6_MA_img, Write_7_VT_img;	
char  Make__5_ID_vid, Make__6_MA_vid, Make__7_VT_vid;		// Flags set in command file to allow ImageJ macro to make videos and then delete the bitmaps used, unless the corresponding Write_... flags are true
char  Wound_Healing_Mode;									// Flag used if only the first cells in the healing front are to be included among valid tracks 
char  Color_Print_Out_Mode = FALSE;							// Flag used for making outlines brighter, while color dots and tracks darker so images look better when printed on paper
int   Top____Crop_Margin;									// Number of pixels to remove from  top   of image (where black stripes may be present) 
int   Bottom_Crop_Margin;									// Number of pixels to remove from bottom of image (where scale bar and time stamp may be present) 
int   Left___Crop_Margin;									// Number of pixels to remove from  left  of image (where black stripes may be present) 
int   Right__Crop_Margin;									// Number of pixels to remove from  right of image (where black stripes may be present) 
int	  Min_Diameter_px;										// Diameter [pixels] of smallest cell to include in analysis
int	  Max_Diameter_px;										// Diameter [pixels] of largest  cell to include in analysis
int   Mid_Diameter_px;										// The midpoint diameter [pixels] between the Min and Max diameter defined above. 
int   Min_Cell_Pixels;										// Smallest number of pixels in a segmented region to keep as an identified cell 
int   Max_Cell_Pixels;										// Largest  number of pixels in a segmented region to keep as an identified cell 
int   Mid_Cell_Pixels;										// Midpoint number of pixels between Min and Max given above 
int   Previous_Image_Shift_X_px, Previous_Image_Shift_Y_px;	// The shift vector which gave best fit for the previous image. Used as starting point for new shift search.
int   First_Image_Number;									// The sequence number (typically 0 or 1) of the first image, given in the file celltraxx_bitmaps2process.txt
int   Last__Image_Number;									// The sequence number of the last image, given in the file celltraxx_bitmaps2process.txt 
int   Image_Nr_Increment;									// The number added to the previous image number to get the next number to process, given in the file celltraxx_bitmaps2process.txt 
int   Tuning_Image_Code;									// 1 if tuning first image in first video, 2 if tuning middle image, 3 if tuning last image
int   Tuning_First_Image;									// The first image number to display when tuning. Also used when determining time span of data in the file _Identified_cells_summary.csv
int   Tuning__Last_Image;									// The  last image number to display when tuning. Also used when determining time span of data in the file _Identified_cells_summary.csv
int   Shortest_Cell_Track;									// The lowest number of levels a matched cell needs to have to be included in the output results and statistics. This removes short lived cells and dust 
int   Gaussian_Filter_Halfwidth;							// The matrix size of the gaussian filter is 2*Gaussian_Filter_Halfwidth+1 
int   Nof_Matched_Cells;									// The currently highest number of matched cells, needed when adding new cells after divisions etc.
int   Nof_Bins;												// Number of bins in output histogram				
int   Track_Smoothing_Iterations;							// Number of times the filtering of track paths will be iterated 	
int   WH_Box_Upper_Limit_px, WH_Box_Lower_Limit_px;			// Upper and lower limit for cells to include during wound healing analysis
int   WH_Nof_Cells_Upper_First, WH_Nof_Cells_Lower_First;	// Number of identified cells in upper and lower part of the first image (level==0), used in wound healing analysis.
int   WH_Nof_Cells_Upper[MAXLEVELS];						// Number of identified cells in upper part of image, used in wound healing analysis. Need to store the number at each level to correctly calculate averages with the right weights for velocity matrix.
int   WH_Nof_Cells_Lower[MAXLEVELS];						// Number of identified cells in lower part of image 
long  Nof_Identified_Cells[MAXLEVELS];						// Array with initial number of identified cells in each image at each level
long  Nof__Cells[MAXLEVELS];								// Array with total number of identified cells in each image at each level, which might be higher then Nof_Identified_Cells[] if the matching interpolation adds a new cell in the previous layer 
long  PointList_x[MAXX*MAXY];								// Used as buffer in the FILL functions
long  PointList_y[MAXX*MAXY];
float TimeStep;												// The time interval [minutes] between consecutive images, used for calculating cell velocity correctly 
float Max_Cell_Velocity;									// Highest allowed cell velocity [µm/min]. Also used as the endpoint of the horizontal axis of the velocity histogram
float Max_Cell_Shift_px;									// Longest allowed distance [pixels] for a cell to move from one image to the next. Used to rule out fast moving contamination and to limit the search for matching cells. 
float Pixel_Size;											// Common square pixel size in all images [µm/pixel]					
float    Dot_Diameter_um;									// Diameter of the tracking color dots 
float    Min_Diameter_um;									// Typical cell diameter [µm] during cell division (when it is smallest) 
float    Max_Diameter_um;									// The largest cell diameter [µm] to include in the analysis. Can be used to exclude unwanted dark regions from images. 
float Max_Image_Shift_um;									// Longest distance [µm] to search when aligning the current image to the previous one when correcting for image misalignment 
int   Cutting_Diameter_px;									// The largest product in pixels of distance and one-sided radius of curvature for cutting a connected region
float Cutting_Diameter_um;									// The diameter [µm] which governs how small cell protrusions will be cut. Zero will omit the cut_concave_region() function.
float Version;												// Must match with version in command file	
float Segmentation_Limit;									// Number of standard deviations as color segmentation interval 
float     Segmentation_Grey_Limit_From_First_Image;	// The grey level in the background (when cells are masked) of the smoothed image which gives a 0.3 % count rate comparded to the most common grey level. Used for automatic segmentation. 
float Previous_Segmentation_Grey_Limit_From_First_Image;	// The previous value written to the command file to be re-used if the same image is read again, to save time and allow tuning of middle and last image to be done using the auto segmentation limit from the first image in the video - as is done when processing videos 
float VT_Image_Contrast;									// Factor used to enhance (>1) or reduce (<1) the contrast in the valid tracks images to be used in videos
float          Gaussian_Radius;								// Radius in µm for gaussian filter
float          Gaussian_Filter[MAXX][MAXY];					// Used in gaussian_smooth() function and filled in gaussian_populate_filter() function 
float  Previous_Expanded_Image[MAXX][MAXY];					// Array for storing the expanded image from the previous level, used in align_images()  
float           Expanded_Image[MAXX][MAXY];					// Array for storing floating point image after mirroring margins for easy smoothing and alignment 
float           Smoothed_Image[MAXX][MAXY];					// Array for storing smoothed image and doing more accurate segmentation with floating point precision 
short          Segmented_Image[MAXX][MAXY];					// Array for storing cell image after segmentation as 0 = background and 1 = cell  
short Previous_Segmented_Image[MAXX][MAXY];					// Array for storing previous cell image after segmentation as 0 = background and 1 = cell  
short             Tuning_Image[MAXX][MAXY];					// Array for storing image data for tuning image (also used as temporary storage) 
short            Matched_Image[MAXX][MAXY];					// Array for storing numbered cells where the numbers are matched to the previous level image
int      WH_Matched_Histogram[     2   ][MAXY];				// Arrays for storing the number of matched cells versus current center coordinate y to allow including only a given number of tracks closest to the image center during Wound Healing (WH) 
int        Diameter_Histogram[     2   ][MAXBINS];			// Array for storing the equivalent cell diameter histograms. Only using three columns in output file: Equivalent diameter bin midpoint [µm], count of all cell candidates and count of all accepted cells 
int        Velocity_Histogram[MAXLEVELS][MAXBINS];			// Array for storing the cell velocity histograms
int             Matched_Cells[MAXLEVELS][MATCHED_CELLS];	// Matrix which summarizes the identification process by showing which cell numbers that were matched in the various levels 
float         Velocity_Matrix[MAXLEVELS][MATCHED_CELLS];	// Array for storing the cell velocities at each level and for each matched cell. 	
float        Image_Shift_um_X[MAXLEVELS];					// The x coordinate of the image shift calculated by align_images() and used for writing this shift to file  
float        Image_Shift_um_Y[MAXLEVELS];					// The y coordinate of the image shift calculated by align_images() and used for writing this shift to file  
float Cell_Shift_Vector_Sum_X[MAXLEVELS];					// The x coordinate of the accumulated shift for all matched cells in a level. Used for calculating the average length of this vector (converted to a velocity) which will be small if there is little image shift, but large if all cells move a bit extra in a certain direction.
float Cell_Shift_Vector_Sum_Y[MAXLEVELS];					// The y coordinate of the accumulated shift for all matched cells in a level. Used for calculating the average length of this vector (converted to a velocity) which will be small if there is little image shift, but large if all cells move a bit extra in a certain direction.
float Cell_Shift_Vector_Sum_X_WH_Upper[MAXLEVELS];			// The x coordinate of the accumulated shift for the upper cells in wound healing mode. Not really used since the movement in the x direction should be close to zero on average. 
float Cell_Shift_Vector_Sum_Y_WH_Upper[MAXLEVELS];			// The y coordinate of the accumulated shift for the upper cells in wound healing mode. Used for calculating the overall cell velocity into the wound.
float Cell_Shift_Vector_Sum_X_WH_Lower[MAXLEVELS];			// The x coordinate of the accumulated shift for the lower cells in wound healing mode. Not really used since the movement in the x direction should be close to zero on average. 
float Cell_Shift_Vector_Sum_Y_WH_Lower[MAXLEVELS];			// The y coordinate of the accumulated shift for the lower cells in wound healing mode. Used for calculating the overall cell velocity into the wound.
float  Matching_Matrix[MAXCELLS+1][MAXCELLS];				// Matrix used when pairing cells in one level with the next level. The extra first index is the row which marks if a cell in the previous level has been matched with a cell in the current level. 
short  Edge_Point_Matrix[MAXX][MAXY];													// Some arrays declared as global and used in cut_concave_region() and radius_of_curvature_for_fitted_circle()
int    Edge_Point_X_px[CUT_MAX_EDGE_POINTS], Edge_Point_Y_px[CUT_MAX_EDGE_POINTS];		// Coordinates of edge points when numbered from the first point at the lower left and then in a clockwise direction
double Edge_Point_Distances[CUT_MAX_EDGE_POINTS][CUT_MAX_EDGE_POINTS];					// The distances between all points on the edge where only half the matrix gest filled due to the symmetry when two points are swapped
double Pi=3.141592653;
char imgext[4]="bmp";
CELL  Cells[MAXLEVELS][MAXCELLS];							// Array holding all identified cells in all levels
FILE *InFile, *SummaryFile, *CellCountFile, *BMPfile;		// Pointers to the input and output text files							
FILE *PosFile, *MatrixFile, *CPYfile, *TuneFile;			// Pointers to the input and output text files							
FILE *CMDWindowFile, *ErrorMessageFile;						// Pointers to the file which contains the table info also written to the command window during a run and the file with error messages to be read by ImageJ macro
FILE *ShiftFile, *TraxFile;									// Pointers to other output text files	
FILE *VelHistFile, *VelMatrFile;							// Pointers to the velocity histogram and velocity matrix output files 	
FILE *DiaHistFile, *GreyHistogramsFile;						// Pointers  to the cell diameter histogram output file	and grey level histogram file
BITMAP Readinimage,   *Rimage;								// The original image and a pointer to it	
BITMAP Bitmapimage,   *Bimage;								// The input  image and a pointer to it		
BITMAP Outputimage,   *Oimage;								// The output image and a pointer to it	
BITMAP Dummyimage,    *Dimage;								// Dummy image used in filtering routine						
BITMAP Scalebarimage, *Simage;								// Predefined image showing the desired scale bar info, matching the pixel size (in nm)						



//******************************  Declaration of functions  **************************************

void read_command_file_header(void);
char read_image_series_info_from_bitmaps2process(void);
void read_scalebar_image(void); 
char read_and_crop_images(int level); 
void align_images(int level); 
void just_write_dummy_shifted_images(int level); 
float average_abs_grey_level_difference(int i, int j);
void smooth_images_and_identify_cells(int level); 
char cut_concave_region(int start_x, int start_y, int level);
double radius_of_curvature_for_fitted_circle(int i);
double mind(double a, double b); 
void write_interactively_tuned_image(void); 
void match_cells_with_previous_level(int level);
void smooth_tracks(int level);
void generate_valid_tracks_images(int level);
void finish_statistics_files(int level); 
void reset_arrays(void);
void draw_scalebar(void); 
void write_error_message(FILE *fp, char message_type, char *message);
int  gaussian_populate_filter(int filterradius);
void gaussian_smooth(int gauss_halfsize);
int  fill_segmented_box(int fillnumber, int withnumber, int startx, int starty, int *x_min, int *x_max, int *y_min, int *y_max);  
int  fill_segmented(    int fillnumber, int withnumber, int startx, int starty);  
int  fill_tuning(       int fillnumber, int withnumber, int startx, int starty);  
void cellnumber2color(     short cellnumber, unsigned char color[3]);
void cellnumber2monocolor( short cellnumber, unsigned char color[3], char colorcode);
char read_bitmap(BITMAP *bimp, char *infilename);
char out_bitmap(BITMAP *bimp, char *outfilename);
void copy_bitmap(BITMAP *frombimp, BITMAP *tobimp);
void convert2color(BITMAP *bimp);
char crop_bitmap(BITMAP *bimp, int top, int bottom, int left, int right);
char paste_grayscale_with_white_transparent(BITMAP *frombimp, BITMAP *tobimp, unsigned char paste_color[3], int startx, int starty);
char translate_color(char *text, unsigned char color[3]);
char splitbmpname(char *textline, BITMAP *bimp);
int	 maxi(int a, int b);
int	 mini(int a, int b);
double absd(double a);

//************************************* Main program  ***********************************

void main()
{
	char message[STRING_LENGTH];
	int image_level;

	read_command_file_header();

	read_scalebar_image(); 

	while( read_image_series_info_from_bitmaps2process() ){

		image_level = 0; 

		while( image_level<MAXLEVELS && read_and_crop_images(image_level) ){

			if( Perform_Shifting  ) align_images(image_level); 
			else just_write_dummy_shifted_images(image_level);

			smooth_images_and_identify_cells(image_level);

			if( Perform_Tuning ) { write_interactively_tuned_image(); fprintf(ErrorMessageFile, "NO ERRORS\n"); exit(0); }

			match_cells_with_previous_level(image_level);

			image_level++;
		}

		if( image_level == MAXLEVELS ) { sprintf(message, "Maximum number of images in video reached (%d).\nPlease increase MAXLEVELS and recomplie celltraxx.c.\n", MAXLEVELS); write_error_message(ErrorMessageFile, ERROR, message); }

		if( Track_Smoothing_Iterations > 0 ) smooth_tracks(image_level);

		finish_statistics_files(image_level);

		if( Write_7_VT_img || Make__7_VT_vid ) generate_valid_tracks_images(image_level);

		reset_arrays(); 		
	}

	fprintf(ErrorMessageFile, "NO ERRORS\n");	// Reporting no errors when aborting normally to allow ImageJ macro to know if celltraxx.exe crashed 
}

//******************************* Definition of functions ********************************

 
// Reads input data from command file and writes most of it back to screen as a confirmation	

void read_command_file_header(void)
{
	char ch, data_drive_letter, infilename[STRING_LENGTH], outfilename[STRING_LENGTH], cpyfilename[STRING_LENGTH], line[STRING_LENGTH], string[STRING_LENGTH], message[STRING_LENGTH];
	int n, r; 
	float vers;	
	
	// Initializing some  of the global variables 
	Version=(float)VERSION;
	Rimage=&Readinimage;
	Bimage=&Bitmapimage;
	Oimage=&Outputimage;
	Dimage=&Dummyimage;
	Simage=&Scalebarimage;

	// Opening the text file for writing error messages to allow the ImageJ macro to display them even after the command window closes upon exit 
	sprintf(outfilename, "C:/celltraxx_system/celltraxx_error_messages.txt");
	if( !(ErrorMessageFile=fopen(outfilename, "w")) ){ printf("Can't create error message file: %s\nPlease check if the file is open in another application.\n", outfilename); exit(0); }

	sprintf(infilename, "C:/celltraxx_system/celltraxx_defaults.txt"); 
	if( (InFile=fopen(infilename, "r"))==NULL ) write_error_message(ErrorMessageFile, ERROR, "Can't open command file celltraxx_defaults.txt.\nPlease check that it is present in the directory C:/celltraxx_system/");	// write_error_message(ErrorMessageFile, FALSE, message);

	// Reading data from command file 
	printf("\n");
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Wound healing mode %s",             &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Wound healing mode'."); }																            else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Wound_Healing_Mode =TRUE; printf("Will analyse images in wound healing mode.        \n"); } else { Wound_Healing_Mode =FALSE; printf("Will not use wound healing mode.              \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Perform flat field correction %s",  &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Perform flat field correction'."); }																else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Perform_FFC        =TRUE; printf("Will work with flat field corrected images.       \n"); } else { Perform_FFC        =FALSE; printf("Will use images without flat field correction.\n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Perform image shift correction %s", &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Perform image shift correction'."); }																else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Perform_Shifting   =TRUE; printf("Will perform image shift correction.              \n"); } else { Perform_Shifting   =FALSE; printf("Will not perform image shift correction.      \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Perform interactive tuning %s",     &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Perform interactive tuning'."); }																	else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Perform_Tuning     =TRUE; printf("Will perform interactive tuning of parameters.    \n"); } else { Perform_Tuning     =FALSE; printf("Will use default parameters from last run.    \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Version %f",                        &vers)<1 || vers!=Version )											{ sprintf(message, "Missing or wrong version number %3.1f. Current version is %3.1f.", vers, Version);																	write_error_message(ErrorMessageFile, ERROR, message); }
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Track smoothing iterations %d",     &Track_Smoothing_Iterations)<1 || Track_Smoothing_Iterations<0 )	{ sprintf(message, "Missing or invalid track smoothing iterations value %d. Should not be negative.", Track_Smoothing_Iterations);										write_error_message(ErrorMessageFile, ERROR, message); }	else { if( Track_Smoothing_Iterations == 0 ) printf("Will not do any smoothing of tracks.\n"); else printf("Will filter smooth tracks using %d iteration(s).\n", Track_Smoothing_Iterations); }
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Results folder drive letter %c",    &ch)<1 || ch<'A' || ch>'Z')											{ sprintf(message, "Missing or invalid drive letter '%c' in command file. Using default drive: C://", ch); data_drive_letter = 'C';										write_error_message(ErrorMessageFile, WARNING, message); }	else { data_drive_letter = ch; printf("Will read image data from folder %c:/celltraxx_data/ \n", data_drive_letter); }  
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Results folder name %s",            &string)<1 )														{ sprintf(message, "Missing or invalid name for results folder. Using default: Results"); sprintf(ResultsFolderName,"%c:/celltraxx_data/Results", data_drive_letter);	write_error_message(ErrorMessageFile, WARNING, message); }	else { sprintf(ResultsFolderName, "%c:/celltraxx_data/%s", data_drive_letter, string); printf("Will write results to sub-folder %s\n", ResultsFolderName); }  
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "First part of folder name %s",      &string)<1 );														// Not used in celltraxx.exe, just added for convenience in ImageJ macro to remember last user setting 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Pixel size [um] %f",	            &Pixel_Size)<1 || Pixel_Size<=0 )									{ sprintf(message, "Missing or invalid pixel size %5.3f.", Pixel_Size);																									write_error_message(ErrorMessageFile, ERROR, message); }	else printf("Pixel size: %7.5f um/pixel.\n", Pixel_Size);
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Gaussian filter radius [um] %f",    &Gaussian_Radius)<1 || Gaussian_Radius<0 )							{ sprintf(message, "Missing or invalid radius for gaussian filter %1.0f pixels.\n", Gaussian_Radius/Pixel_Size );														write_error_message(ErrorMessageFile, ERROR, message); }	else printf("Smoothing images with Gaussian filter of radius  %4.1f um = %d pixels.\n", Gaussian_Radius, r=(int)(Gaussian_Radius/Pixel_Size+0.5) ); 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Smallest cell diameter [um] %f",    &Min_Diameter_um)<1 || Min_Diameter_um/Pixel_Size<1 )				{ sprintf(message, "Missing or invalid smallest cell diameter %4.2f pixels. Should be above 1 pixels.", Min_Diameter_um/Pixel_Size );									write_error_message(ErrorMessageFile, ERROR, message); }	else printf("Analysing only cells with diameters larger  than %4.1f um.\n", Min_Diameter_um ); 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Largest  cell diameter [um] %f",    &Max_Diameter_um)<1 || Max_Diameter_um<=Min_Diameter_um )			{ sprintf(message, "Largest cell diameter (%3.1f pixels) < smallest cell diameter (%3.1f pixels).", Max_Diameter_um/Pixel_Size, Min_Diameter_um/Pixel_Size );			write_error_message(ErrorMessageFile, ERROR, message); }	else printf("Analysing only cells with diameters smaller than %4.1f um.\n", Max_Diameter_um ); 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Cutting  cell diameter [um] %f",    &Cutting_Diameter_um)<1 || Cutting_Diameter_um<0 )	   				{ sprintf(message, "Missing or invalid cutting cell diamter %4.2f.", Cutting_Diameter_um);																			    write_error_message(ErrorMessageFile, ERROR, message); }	else if( Cutting_Diameter_um == 0 ) printf("Will not cut cell candidate regions at narrow, concave points since 'Cutting cell diameter' = 0.\n"); else printf("Will cut cell candidate regions at narrow, concave points if the cutting distance is < %3.1f um.\n", Cutting_Diameter_um);
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Top    crop margin [pixels] %d",    &Top____Crop_Margin)<1 || Top____Crop_Margin < 0 )					{ sprintf(message, "Top    crop margin %d. Must be zero or larger.", Top____Crop_Margin);																				write_error_message(ErrorMessageFile, ERROR, message); }	else printf("Will crop off %d pixels from  top   of all images.\n", Top____Crop_Margin);
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Bottom crop margin [pixels] %d",    &Bottom_Crop_Margin)<1 || Bottom_Crop_Margin < 0 )					{ sprintf(message, "Bottom crop margin %d. Must be zero or larger.", Bottom_Crop_Margin);																				write_error_message(ErrorMessageFile, ERROR, message); }	else printf("Will crop off %d pixels from bottom of all images.\n", Bottom_Crop_Margin);
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Left   crop margin [pixels] %d",    &Left___Crop_Margin)<1 || Left___Crop_Margin < 0 )					{ sprintf(message, "Left   crop margin %d. Must be zero or larger.", Left___Crop_Margin);																				write_error_message(ErrorMessageFile, ERROR, message); }	else printf("Will crop off %d pixels from  left  of all images.\n", Left___Crop_Margin);
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Right  crop margin [pixels] %d",    &Right__Crop_Margin)<1 || Right__Crop_Margin < 0 )					{ sprintf(message, "Right  crop margin %d. Must be zero or larger.", Right__Crop_Margin);																				write_error_message(ErrorMessageFile, ERROR, message); }	else printf("Will crop off %d pixels from  right of all images.\n", Right__Crop_Margin);
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Time between images      [min] %f", &TimeStep)<1 || TimeStep<=0 )										{ sprintf(message, "Missing or invalid time step between images %1.1f minutes.", TimeStep );																			write_error_message(ErrorMessageFile, ERROR, message); }	else printf("The time step between images is %3.1f min.\n", TimeStep ); 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Highest cell velocity [um/min] %f", &Max_Cell_Velocity)<1 || Max_Cell_Velocity<0 )						{ sprintf(message, "Missing or invalid value for highest cell velocity %1.1f um/min.", Max_Cell_Velocity );																write_error_message(ErrorMessageFile, ERROR, message); }	else printf("The highest allowed cell velocity is %1.1f um/min, corresponding to a movement of %d pixels between consecutive images.\n", Max_Cell_Velocity, (int)(Max_Cell_Velocity*TimeStep/Pixel_Size+0.5) ); 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Shortest cell track   [images] %d", &Shortest_Cell_Track)<1 || Shortest_Cell_Track<3 )					{ sprintf(message, "Missing or invalid shortest cell track %d. Should cover at least 3 images.", Shortest_Cell_Track);													write_error_message(ErrorMessageFile, ERROR, message); }	else printf("Including in the results files only cells that were tracked over at least %d steps (%d consecutive images).\n", Shortest_Cell_Track, Shortest_Cell_Track+1);
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "First image number %d",				&Tuning_First_Image)<1 || Tuning_First_Image<0 )					{ sprintf(message, "Missing or invalid first image number (%d). Can not be negative.\n", Tuning_First_Image);															write_error_message(ErrorMessageFile, ERROR, message); }
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Last  image number %d",				&Tuning__Last_Image)<1 || Tuning__Last_Image<Tuning_First_Image )	{ sprintf(message, "Missing or invalid  last image number (%d). Must be larger than first image.\n", Tuning__Last_Image);												write_error_message(ErrorMessageFile, ERROR, message); }	else printf("Will count cells from image number %d to %d, coresponding to times %3.1f to %3.1f minutes.\n", Tuning_First_Image, Tuning__Last_Image, Tuning_First_Image*TimeStep, Tuning__Last_Image*TimeStep ); 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Image number increment %d",         &Image_Nr_Increment)<1 || Image_Nr_Increment<1 )					{ sprintf(message, "Missing or invalid image number increment %d. Should be a positive number.\n", Image_Nr_Increment);													write_error_message(ErrorMessageFile, ERROR, message); }	else { if( Image_Nr_Increment == 1 ) printf("Will analyse consecutive images.\n"); else printf("Will analyze every %d images.\n", Image_Nr_Increment); }
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Make identified cell videos %s",    &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Make identified cell videos'."); }																	else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Make__5_ID_vid     =TRUE; printf("Will tell ImageJ to make identified cell videos.  \n"); } else { Make__5_ID_vid     =FALSE; printf("Will not make identified cell videos.         \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Make matched cell videos %s",       &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Make matched cell videos'."); }																	else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Make__6_MA_vid     =TRUE; printf("Will write matched cell images for checking.      \n"); } else { Make__6_MA_vid     =FALSE; printf("Will not make matched cell videos.            \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Make valid track videos %s",		&string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Make valid track videos'."); }																		else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Make__7_VT_vid     =TRUE; printf("Will write images with valid tracks for video.    \n"); } else { Make__7_VT_vid     =FALSE; printf("Will not make valid track videos.             \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Tracking dot diameter [um] %f",     &Dot_Diameter_um)<1 || Dot_Diameter_um<0 )							{ sprintf(message, "Missing or invalid tracking dot diameter %4.2f um. Should a positive number.", Dot_Diameter_um );													write_error_message(ErrorMessageFile, ERROR, message); }	else if( Dot_Diameter_um == 0 ) {  printf("Will not draw color dots at center of mass for matched and tracked cells.\n"); Draw_CM=FALSE; } else { printf("Will draw color dots with diameter %3.1f um at center of mass for matched and tracked cells.\n", Dot_Diameter_um ); Draw_CM=TRUE; }
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Valid track image contrast %f",     &VT_Image_Contrast)<1 || VT_Image_Contrast<=0 )						{ sprintf(message, "Missing or invalid image contrast factor %4.2f. Should be a positive number.", VT_Image_Contrast);													write_error_message(ErrorMessageFile, ERROR, message); }    else printf("Valid track image contrast factor: %3.1f (where > 1.0 gives contrast enhancement).\n", VT_Image_Contrast);	
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Scale bar color %s",                &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Scale bar color'."); }																                else { strcpy(Scalebar_Color, string); printf("Will draw scale bar with %s color in images if scale bar image exists.\n", Scalebar_Color);  } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Draw cell outline %s",              &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Draw cell outline'."); }																            else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Draw_Outline       =TRUE; printf("Will draw colored outline in matched images.      \n"); } else { Draw_Outline       =FALSE; printf("Will not mark outline in matched images.      \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Draw cell track line %s",           &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Draw cell track line'."); }																		else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Draw_Track         =TRUE; printf("Will draw colored track of cell movements.        \n"); } else { Draw_Track         =FALSE; printf("Will not mark cell track in output image.     \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Write mirror margin images %s",     &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Write mirror margin images'."); }																	else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Write_0_MM_img     =TRUE; printf("Will write mirror margin images for checking.     \n"); } else { Write_0_MM_img     =FALSE; printf("Will not write mirror margin images.          \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Write shifted images %s",			&string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Write shifted images'."); }																		else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Write_1_SH_img     =TRUE; printf("Will write shifted images for checking.           \n"); } else { Write_1_SH_img     =FALSE; printf("Will not write shifted images.                \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Write gaussian smoothed images %s", &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Write gaussian smoothed images'."); }																else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Write_2_SM_img     =TRUE; printf("Will write smoothed images for checking.          \n"); } else { Write_2_SM_img     =FALSE; printf("Will not write smoothed images.               \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Write segmented cell images %s",    &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Write segmented cell images'."); }																	else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Write_3_SE_img     =TRUE; printf("Will write segmented cell images for checking.    \n"); } else { Write_3_SE_img     =FALSE; printf("Will not write segmented images.              \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Write cut cell images %s",          &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Write cut cell images'."); }																	    else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Write_4_CU_img     =TRUE; printf("Will write cut cell images for checking.          \n"); } else { Write_4_CU_img     =FALSE; printf("Will not write cut cell images.               \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Write identified cell images %s",   &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Write identified cell images'."); }																else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Write_5_ID_img     =TRUE; printf("Will write identified cell images for checking.   \n"); } else { Write_5_ID_img     =FALSE; printf("Will not write identified cell images.        \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Write matched cell images %s",      &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Write matched cell images'."); }																	else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Write_6_MA_img     =TRUE; printf("Will write matched cell images for checking.      \n"); } else { Write_6_MA_img     =FALSE; printf("Will not write matched cell images.           \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Write valid track images %s",       &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Write valid track images'."); }																	else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Write_7_VT_img     =TRUE; printf("Will write images with valid tracks for checking. \n"); } else { Write_7_VT_img     =FALSE; printf("Will not write valid track images.            \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Keep cells from previous image %s", &string)<1 )														{ write_error_message(ErrorMessageFile, ERROR, "Missing command for 'Keep cells from previous image'."); }																else if( !strcmp(string,"yes") || !strcmp(string,"Yes") || !strcmp(string,"YES") ) { Keep_Previous_Cells=TRUE; printf("Will use previous cells to find new cell regions. \n"); } else { Keep_Previous_Cells=FALSE; printf("Will not use previous cells in next image.    \n"); } 
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Segmentation limit [SDs] %f",       &Segmentation_Limit)<1  || Segmentation_Limit==0 )					{ sprintf(message, "Missing or invalid segmentation limit %4.2f. Can not be zero.", Segmentation_Limit);																write_error_message(ErrorMessageFile, ERROR, message); }	if( Segmentation_Limit < 0 ) { Use_SD_Based_Segmentation_Grey_Limit = TRUE; Segmentation_Limit = -Segmentation_Limit; printf("Will use as segmentation limit the smoothed first image mean grey level minus %4.2f * the grey level standard deviation.\n", Segmentation_Limit); } else { Use_SD_Based_Segmentation_Grey_Limit = FALSE; printf("Will automatically find the segmentation limit from the background grey level of the first image\n"); }
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "# bins in histogram %d",            &Nof_Bins)<1 || Nof_Bins>=MAXBINS ||  Nof_Bins<2 ) 					{ sprintf(message, "Missing or invalid number of bins %d. Should be between 2 and %d.", Nof_Bins, MAXBINS-1);															write_error_message(ErrorMessageFile, ERROR, message); }	else printf("Number of bins in histogram:  %d\n", Nof_Bins);
	if( (fgets(line, STRING_LENGTH, InFile)==NULL) || sscanf(line, "Tuning image code %d",              &Tuning_Image_Code)<1 || Tuning_Image_Code<0 || Tuning_Image_Code>3){ sprintf(message, "Missing or invalid tuning image code %d. Should be 0, 1, 2 or 3.\n", Tuning_Image_Code);															write_error_message(ErrorMessageFile, ERROR, message); }	else { if( Tuning_Image_Code  == 1 ) printf("Will display first image during tuning.\n"); else if( Tuning_Image_Code == 2 ) printf("Will display the middle image during tuning.\n"); else printf("Will display the last image during tuning.\n"); }

	sprintf( InputFolderName, "%c:/celltraxx_data", data_drive_letter);					// Reading bitmap images from this folder 
	if( Perform_Tuning ) {
		sprintf(      ResultsFolderName, "C:/celltraxx_system/tuning_temp");					// Writing tuning results data to separate, temporary folder  
		if(      _mkdir(ResultsFolderName) == 0 ) printf("The folder %s was sucessfully created!\n", ResultsFolderName); else printf("Writing temporary results to folder %s\n", ResultsFolderName);
	} else { if( _mkdir(ResultsFolderName) == 0 ) printf("The folder %s was sucessfully created!\n", ResultsFolderName); else printf("Warning: The folder %s already exists. Old data will be overwritten.\n", ResultsFolderName); 
		Tuning_Image_Code = 0;		// Needed to avoid analysing the same images repeatedly if the user clicked [Back] in ImageJ and then unchecked "Perform interactive tuning" before clicking [OK] to run celltraxx.exe directly
	}

	Cutting_Diameter_px = (int)( Cutting_Diameter_um / Pixel_Size ); 

	Max_Cell_Shift_px = Max_Cell_Velocity * TimeStep * Image_Nr_Increment / Pixel_Size;			// Units: pixels = µm/min * min * # timesteps to skip between images / (µm/pixel) 

	Gaussian_Filter_Halfwidth = gaussian_populate_filter(r);					// r is calculated at the end of the input line for "Gaussian filter radius" above. In ImageJ, the Process > Filters > Gaussian Blur filter radius "Sigma (Radius)" is 3.03*Pixel_size times SMALLER than r = Gaussian_Filter_Halfwidth + 1 = (int)(Gaussian_Radius/Pixel_Size+0.5), so  Gaussian_Radius = 3.03 * Pixel_Size * Sigma(Radius) in ImageJ 

	// Reading data from tuning info file (and later writing updated info back to the same file as a way to remember the segmentation limit if the same image is analysed again and if tuning the middle and last image of the same video).
	sprintf(infilename, "C:/celltraxx_system/celltraxx_tuning_info.txt"); 
	if( (TuneFile=fopen(infilename, "r"))==NULL ){ sprintf(message, "Can't open file celltraxx_tuning_info.txt for reading.\nPlease check that it is present in the directory C:/celltraxx_system/"); write_error_message(ErrorMessageFile, ERROR, message); }
	if( (fgets(line, STRING_LENGTH, TuneFile)==NULL) || sscanf(line, "First image grey level limit %f",   &Previous_Segmentation_Grey_Limit_From_First_Image)<1 || Previous_Segmentation_Grey_Limit_From_First_Image<0 )	{ sprintf(message, "Missing or invalid grey level limit %5.3f. Should be a positive number.\nPlease check the value in the file\n%s.\n", Previous_Segmentation_Grey_Limit_From_First_Image, infilename); write_error_message(ErrorMessageFile, ERROR, message);}	else printf("Previous grey level limit: %5.3f\n", Previous_Segmentation_Grey_Limit_From_First_Image);	
	if( (fgets(line, STRING_LENGTH, TuneFile)==NULL) || sscanf(line, "First image path and name    %s",   &PreviousFirstImageName)<1 )																						{ sprintf(message, "Missing name for first image. Please check the name in the file\n%s.\n", infilename); write_error_message(ErrorMessageFile, ERROR, message);}																									else printf("Path and name of first image in previous run:\n%s\n", PreviousFirstImageName);   
	fclose(TuneFile); 

	// Opening common results file with summary data for all videos in celltraxx_bitmaps2process.txt and writing file header 
	r = strlen(infilename); 
	while( --r>3 && infilename[r]!='.' ) infilename[r]='\0'; 
	infilename[r]='\0';													// Removing extension and . in command file file name
	sprintf(outfilename, "%s/_Overall_summary.csv", ResultsFolderName);
	if( !(SummaryFile=fopen(outfilename, "w")) ){ sprintf(message, "Can't create overall summary file \n   %s \nPlease check if the file is open in another application.", outfilename); write_error_message(ErrorMessageFile, ERROR, message); }
	fprintf(SummaryFile, ",Results based on:,,,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Tracks,Time steps,Time steps,Time steps,Time steps,All cell movements,All cell movements,All cell movements,All cell movements,Histogram,Histogram,Histogram,Histogram\n");
	fprintf(SummaryFile, ",,,,Total number,Velocity,Velocity,Velocity,Directness,Directness,Directness,FMI X,FMI X,FMI X,FMI Y,FMI Y,FMI Y,Accumulated distance,Accumulated distance,Accumulated distance,Euclidian distance,Euclidian distance,Euclidian distance,Euclidian velocity,Euclidian velocity,Euclidian velocity,Total number,Velocity,Velocity,Velocity,Total number,Velocity,Velocity,Velocity,Velocity,Velocity,Velocity,Velocity\n");	
	fprintf(SummaryFile, "First image,,Last image,,n,Mean,SD,SEM,Mean,SD,SEM,Mean,SD,SEM,Mean,SD,SEM,Mean,SD,SEM,Mean,SD,SEM,Mean,SD,SEM,n,Mean,SD,SEM,n,Mean,SD,SEM,Mean,Mode,Median (Q2)           (50th percentile),Upper quartile (Q3) (75th percentile),Video number,Well coordinates\n");
	fprintf(SummaryFile,"Name,# cells,Name,# cells,[ ],[um/min],[um/min],[um/min],[ ],[ ],[ ],[ ],[ ],[ ],[ ],[ ],[ ],[um],[um],[um],[um],[um],[um],[um/min],[um/min],[um/min],[ ],[um/min],[um/min],[um/min],[ ],[um/min],[um/min],[um/min],[um/min],[um/min],[um/min],[um/min],,Letter,Number\n");

	// Opening the file for writing the table data also written to the command window 
	sprintf(outfilename, "%s/_Log_window_summary.csv", ResultsFolderName);
	if( !(CMDWindowFile=fopen(outfilename, "w")) ){ sprintf(message, "Can't create log window summary file \n%s \nPlease check if the file is open in another application.", outfilename); write_error_message(ErrorMessageFile, ERROR, message); }
	fprintf(CMDWindowFile, "Image name,Shift X,Shift Y,Grey,Grey,SD based,Applied,# cells,# cells\n");	// Writing header only once at the top two lines in the .csv file 
	fprintf(CMDWindowFile, ",[um],[um],mean,SD,grey limit,grey limit,in image,matched\n");   

	// Opening the file for writing the cell count data to a separate file  
	sprintf(outfilename, "%s/_Identified_cells_summary.csv", ResultsFolderName);
	if( !(CellCountFile=fopen(outfilename, "w")) ){ sprintf(message, "Can't create cell count summary file \n%s \nPlease check if the file is open in another application.", outfilename); write_error_message(ErrorMessageFile, ERROR, message); }
	fprintf(CellCountFile,                                           ",,,,Image number,");		for( n=0; n<=Tuning__Last_Image; n++ ) fprintf(CellCountFile, "%d,",	n*Image_Nr_Increment              );		fprintf(CellCountFile, "\n"); 
	fprintf(CellCountFile,                                           ",,,,Time [min],"  );		for( n=0; n<=Tuning__Last_Image; n++ ) fprintf(CellCountFile, "%3.1f,", n*Image_Nr_Increment*TimeStep     );		fprintf(CellCountFile, "\n"); 
	fprintf(CellCountFile,                                           ",,,,Time [h],"    );		for( n=0; n<=Tuning__Last_Image; n++ ) fprintf(CellCountFile, "%5.3f,", n*Image_Nr_Increment*TimeStep/60  ); 		fprintf(CellCountFile, "\n"); 
	fprintf(CellCountFile, "Video name,VID number,Well letter,Well number,Time [days]," );		for( n=0; n<=Tuning__Last_Image; n++ ) fprintf(CellCountFile, "%8.6f,", n*Image_Nr_Increment*TimeStep/3600); 		fprintf(CellCountFile, "\n"); 

	// Opening the file for writing the grey level histogram data 
	sprintf(outfilename, "%s/_Grey_level_curves_for_first_images.csv", ResultsFolderName);
	if( !(GreyHistogramsFile=fopen(outfilename, "w")) ){ sprintf(message, "Can't create the file with grey level curves for first images \n   %s\nPlease check if the file is open in another application.", outfilename); write_error_message(ErrorMessageFile, ERROR, message); }
	fprintf(GreyHistogramsFile, "Image name,Bin number,Grey value,Image pixel count,Background pixel count,Log10(Normalized pixel count),Fitted parabola\n");	// Writing header only once at the top two lines in the .csv file 
	
	// Copying the current celltraxx_bitmaps2process.txt file into the given results directory 
	if( (BMPfile=fopen("C:/celltraxx_system/celltraxx_bitmaps2process.txt", "r"))==NULL ){ sprintf(message, "Can not find the file celltraxx_bitmaps2process.txt in this folder.\nPlease check if the file exists.\n");  write_error_message(ErrorMessageFile, ERROR, message); }
	sprintf(cpyfilename, "%s/celltraxx_bitmaps2process.txt", ResultsFolderName);
	if( (CPYfile=fopen(cpyfilename, "w"))==NULL ){ sprintf(message, "Can not copy the file celltraxx_bitmaps2process.txt to the folder \n   %s \nPlease check if the file is open in another application.\n", ResultsFolderName);  write_error_message(ErrorMessageFile, ERROR, message); }
	ch = fgetc(BMPfile);
	while(ch != EOF) {
		fputc(ch, CPYfile);
		ch = fgetc(BMPfile);
	}
	printf("Copied the file celltraxx_bitmaps2process.txt to the current results folder: %s\n", ResultsFolderName);
	fclose(CPYfile);
	rewind(BMPfile);	// After copying file content, going back to top and checking that header is present.
	if( (fgets(line, STRING_LENGTH, BMPfile)==NULL) || sscanf(line, "%s", &BMPfilename)<1 || strcmp(BMPfilename, "Bitmap") ) { sprintf(message, "Missing header line in file celltraxx_bitmaps2process.txt\nJust found: %s", line);  write_error_message(ErrorMessageFile, ERROR, message); } 

	// Copying the current command file into the given results directory for future reference
	sprintf(cpyfilename, "%s/celltraxx_settings.txt", ResultsFolderName);
	if( (CPYfile=fopen(cpyfilename, "w"))==NULL ){ sprintf(message, "Can not copy the file celltraxx_defaults.txt to the folder \n   %s \nwith new name celltraxx_settings.txt \nPlease check if the file is open in another application.", ResultsFolderName);  write_error_message(ErrorMessageFile, ERROR, message); }
	rewind(InFile);
	ch = fgetc(InFile);
	while(ch != EOF) {
		fputc(ch, CPYfile);
		ch = fgetc(InFile);
	}
	printf("Copied the currently used settings to the results folder with name: \n%s\n", cpyfilename);
	fclose(CPYfile);
	fclose(InFile);
}


// Reading the next first file name and last image number from the bitmap file. Returning TRUE if successful, FALSE if no more image names found. 

char read_image_series_info_from_bitmaps2process(void)
{
	char outfilename[STRING_LENGTH], line[STRING_LENGTH], message[STRING_LENGTH];

	if( (fgets(line, STRING_LENGTH, BMPfile)==NULL) || sscanf(line, "%s %d %d", &BMPfilename, &First_Image_Number, &Last__Image_Number)<3 ) { printf("\nFound no (more) valid image series names in the file celltraxx_bitmaps2process.txt. Last line is:\n%s\n", line); return(FALSE); } else printf("\nFirst image to read is: %s%04d.bmp \n", BMPfilename, First_Image_Number);

	if( First_Image_Number+Shortest_Cell_Track > Last__Image_Number+1) { sprintf(message, "The last image number (%04d) in the file celltraxx_bitmaps2process.txt is too small. \nIt must be at least %d to give any data, since the first image number is %d, and the shortest cell track is %d images.\nPlease update the settings.\n", Last__Image_Number, First_Image_Number+Shortest_Cell_Track, First_Image_Number, Shortest_Cell_Track); write_error_message(ErrorMessageFile, ERROR, message); } else printf("Last  image to read is: %s%04d.bmp \n", BMPfilename, Last__Image_Number);
	if( First_Image_Number+Image_Nr_Increment  > Last__Image_Number  ) { sprintf(message, "The last image number (%04d) in the file celltraxx_bitmaps2process.txt is too small when the first image number is %04d and the step between images is %04d.\nPlease update the settings.\n", Last__Image_Number, First_Image_Number, Image_Nr_Increment); write_error_message(ErrorMessageFile, ERROR, message); }                                                                                          // else printf("Will analyse every %d image(s).\n", Image_Nr_Increment);

	// Opening the file containing the cell variable values 
	sprintf(outfilename, "%s/%s_Position-data.csv", ResultsFolderName, BMPfilename);
	if( !(PosFile=fopen(outfilename, "w")) ){ sprintf(message, "Can't create output data file: %s\nPlease check if the file is open.\n", outfilename); write_error_message(ErrorMessageFile, ERROR, message); }

	// Opening the file containing the histogram data 
	sprintf(outfilename, "%s/%s_Matched-cells.csv", ResultsFolderName, BMPfilename);
	if( !(MatrixFile=fopen(outfilename, "w")) ){ sprintf(message, "Can't create output data file: %s\nPlease check if the file is open.\n", outfilename); write_error_message(ErrorMessageFile, ERROR, message); }

	// Opening the file containing the track specific data 
	sprintf(outfilename, "%s/%s_Track-data.csv", ResultsFolderName, BMPfilename);
	if( !(TraxFile=fopen(outfilename, "w")) ){ sprintf(message, "Can't create output data file: %s\nPlease check if the file is open.\n", outfilename); write_error_message(ErrorMessageFile, ERROR, message); }

	// Opening the file containing the diameter histogram for all matched cells 
	sprintf(outfilename, "%s/%s_Diameter-histogram.csv", ResultsFolderName, BMPfilename);
	if( !(DiaHistFile=fopen(outfilename, "w")) ){ sprintf(message, "Can't create output data file: %s\nPlease check if the file is open.\n", outfilename); write_error_message(ErrorMessageFile, ERROR, message); }

	// Opening the file containing the velocity histogram for all steps of all matched cells 
	sprintf(outfilename, "%s/%s_Velocity-histogram.csv", ResultsFolderName, BMPfilename);
	if( !(VelHistFile=fopen(outfilename, "w")) ){ sprintf(message, "Can't create output data file: %s\nPlease check if the file is open.\n", outfilename); write_error_message(ErrorMessageFile, ERROR, message); }

	// Opening the file containing the velocity histogram for all steps of all matched cells 
	sprintf(outfilename, "%s/%s_Velocity-matrix.csv", ResultsFolderName, BMPfilename);
	if( !(VelMatrFile=fopen(outfilename, "w")) ){ sprintf(message, "Can't create output data file: %s\nPlease check if the file is open.\n", outfilename); write_error_message(ErrorMessageFile, ERROR, message); }

	// Opening the file containing the image shift matrices for checking 
	sprintf(outfilename, "%s/%s_Image-shift-matrices.csv", ResultsFolderName, BMPfilename);
	if( !(ShiftFile=fopen(outfilename, "w")) ){ sprintf(message, "Can't create output data file: %s\nPlease check if the file is open.\n", outfilename); write_error_message(ErrorMessageFile, ERROR, message); }

	return(TRUE); 
}



// Opens a pre made image with text "100 µm" if it exists. If so, the image will be pasted with transparity into the lower right corner of the final Valid Tracks images, above the 100 µm long scale bar

void read_scalebar_image(void)
{
	char imagename[STRING_LENGTH];
	unsigned char col[3];

	Add_Scalebar = Add_100um_image = FALSE;

	sprintf(imagename, "C:/celltraxx_system/celltraxx_100um.bmp");									// Reading the pre defined image with text "100 µm"
	if( translate_color(Scalebar_Color, col) ) {		
		printf("Will draw a 100 um long, %s scale bar on all images.\n", Scalebar_Color); 
		Add_Scalebar = TRUE; 
		if( read_bitmap(Simage, imagename) ) {		
			printf("Will add the text '100 um' in %s above the scale bar.\n", Scalebar_Color); 
			Add_100um_image = TRUE; 
		} 
	} else printf("No scale bars will be added to images.\n"); 
}


// Reads the next image from file, applies the crop margins and shifts the image data diagonally in the matri and mirrors the pixels at the crop margins to allow easy filtering below

char read_and_crop_images(int level) 
{
	char imagename[STRING_LENGTH], dummy_string[3];
	int  i, n, x, y, namelength, top_margin, right_margin; 
	double greyvalue;

	if( First_Image_Number+Image_Nr_Increment*level > Last__Image_Number ) { printf("Finished processing the last image: %s%04d.bmp.\n", BMPfilename, First_Image_Number+Image_Nr_Increment*(level-1)); return(FALSE); } 		// Reading images until last given image number is passed
	if(      Tuning_Image_Code == 1 ) n = First_Image_Number;
	else if( Tuning_Image_Code == 2 ) n = (int)( 0.5*(First_Image_Number+Last__Image_Number) );		// Displaying the middle image during tuning 
	else if( Tuning_Image_Code == 3 ) n = Last__Image_Number;
	else                              n = First_Image_Number+Image_Nr_Increment*level; 
	sprintf(imagename, "%s/%s%04d.bmp", InputFolderName, BMPfilename, n);
	dummy_string[0]='\0';	// Just an empty string to use below 

	if( read_bitmap(Rimage, imagename) ) {					// ...or no more images in the series are found 

		namelength = strlen(Rimage->filename);

		if( level%30 == 0 ){								// Writing the header to screen every 20 levels so the info is always visible 
			printf("\n");   
			printf("%*s", namelength-10, dummy_string); printf("              Shift X   Shift Y   Grey    Grey     SD based       Applied     # cells    # cells \n");   
			printf("%*s", namelength-10, dummy_string); printf("Image name     [um]      [um]     mean     SD     grey limit     grey limit   in image   matched \n");   
		}

		// Making a working copy in Bimage for cropping etc. 
		copy_bitmap(Rimage, Bimage); 

/*		// Simple smooting filter to make images with high contrast regions look more similar if some are shifted much and thus pixels are mixed more in the fractional shifting above. 
		for( y=1; y<Bimage->ysize-1; y++ )
			for( x=1; x<Bimage->xsize-1; x++ ) {
				greyvalue = 0.4*Rimage->image[x][y][0] + 0.15*Rimage->image[x+1][y][0] + 0.15*Rimage->image[x][y+1][0] + 0.15*Rimage->image[x-1][y][0] + 0.15*Rimage->image[x][y-1][0]; 
				for( i=0; i<3; i++) Bimage->image[x][y][i] = (unsigned char)(greyvalue+0.5);  
			}	//*/

		// Cropping images to remove scale bar, time stamp + possibly black stripes on the top and sides 
		crop_bitmap(Bimage, Top____Crop_Margin, Bottom_Crop_Margin, Left___Crop_Margin, Right__Crop_Margin);   

		// Copying image (which now has a new size in x and y) to shifted position in float matrix Expanded_Image[][]
		for( y=0; y<Bimage->ysize; y++ )
			for( x=0; x<Bimage->xsize; x++ )
				Expanded_Image[SHIFTXY+x][SHIFTXY+y] = Bimage->image[x][y][0];				// Using only one channel since original image is 8 bit grey level image. Shifting origin to allow easy filtering and alignment later.

		if( Wound_Healing_Mode && level == 0 ) {											// Setting widest possible wound healing borders for first image, to be narrowed gradually when cell positions are available 						
			WH_Box_Lower_Limit_px = 0;
			WH_Box_Upper_Limit_px = Bimage->ysize-1;
		}

/*		// White border on original image 
		for( y=0; y<Bimage->ysize; y++ )
			Expanded_Image[SHIFTXY][SHIFTXY+y] = Expanded_Image[SHIFTXY+Bimage->xsize-1][SHIFTXY+y]=255;				
		for( x=0; x<Bimage->xsize; x++ )
			Expanded_Image[SHIFTXY+x][SHIFTXY] = Expanded_Image[SHIFTXY+x][SHIFTXY+Bimage->ysize-1]=255;	*/	

		// Mirroring lower part of image to filter/shift margin
		for( y=1; y<=SHIFTXY; y++ )
			for( x=0; x<Bimage->xsize; x++ )
				Expanded_Image[SHIFTXY+x][SHIFTXY-y] = Expanded_Image[SHIFTXY+x][SHIFTXY+y];					

		// Mirroring upper part of image to filter/shift margin
		top_margin = MAXY - Bimage->ysize - SHIFTXY; 
		for( y=1; y<=top_margin; y++ )
			for( x=0; x<Bimage->xsize; x++ )
				Expanded_Image[SHIFTXY+x][SHIFTXY+Bimage->ysize-1+y] = Expanded_Image[SHIFTXY+x][SHIFTXY+Bimage->ysize-1-y];			// Have to set the levels and indeces like this to get the right mirroring about the outer pixel border			

		// Filling all of left margin with mirror image 
		for( y=0; y<MAXY; y++ )
			for( x=0; x<SHIFTXY; x++ )
				Expanded_Image[x][y] = Expanded_Image[2*SHIFTXY-x][y];					

		// Filling all of right margin with mirror image 
		right_margin = MAXX - Bimage->xsize - SHIFTXY; 
		for( y=0; y<MAXY; y++ )
			for( x=0; x<right_margin; x++ )
				Expanded_Image[SHIFTXY+Bimage->xsize+x][y] = Expanded_Image[SHIFTXY+Bimage->xsize-2-x][y];					

		// Writing EXPANDED image to file for checking 
		if( Write_0_MM_img ) {
			copy_bitmap(Bimage, Oimage);		// Just to get header info correct before expanding size 
			Oimage->xsize = (long)(Bimage->xsize + 2*SHIFTXY); if( Oimage->xsize > MAXX ) Oimage->xsize = MAXX;	
			Oimage->ysize = (long)(Bimage->ysize + 2*SHIFTXY); if( Oimage->ysize > MAXY ) Oimage->ysize = MAXY;	
			for( y=0; y<Oimage->ysize; y++ )
				for( x=0; x<Oimage->xsize; x++ )
					for( i=0; i<3; i++ ) Oimage->image[x][y][i]=(unsigned char)(Expanded_Image[x][y]);		// Have to use all three channels to fully overwrite copied Bimage 
			sprintf(imagename, "%s/%s_00_%04d_Mirror-margins.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*level);
			out_bitmap(Oimage, imagename);    //*/
		}
		
//		printf("%32s ", Rimage->filename);
		printf("%*s   ", namelength,    Rimage->filename);
		fprintf(CMDWindowFile, "%s,", Rimage->filename);
		
		// Drawing a dotted lines at the border of the cropped regions
		if( Bottom_Crop_Margin > 0 )	for( x=0; x<Bimage->xsize+1; x+=2 ) Rimage->image[Left___Crop_Margin-1+x            ][Bottom_Crop_Margin-1              ][0]=Rimage->image[Left___Crop_Margin-1+x            ][Bottom_Crop_Margin-1              ][1]=Rimage->image[Left___Crop_Margin-1+x            ][Bottom_Crop_Margin-1              ][2]=0;
		if( Top____Crop_Margin > 0 )	for( x=0; x<Bimage->xsize+1; x+=2 ) Rimage->image[Left___Crop_Margin-1+x            ][Bottom_Crop_Margin  +Bimage->ysize][0]=Rimage->image[Left___Crop_Margin-1+x            ][Bottom_Crop_Margin  +Bimage->ysize][1]=Rimage->image[Left___Crop_Margin-1+x            ][Bottom_Crop_Margin  +Bimage->ysize][2]=0;		
		if( Left___Crop_Margin > 0 )	for( y=0; y<Bimage->ysize+1; y+=2 ) Rimage->image[Left___Crop_Margin-1              ][Bottom_Crop_Margin-1+y            ][0]=Rimage->image[Left___Crop_Margin-1              ][Bottom_Crop_Margin-1+y            ][1]=Rimage->image[Left___Crop_Margin-1              ][Bottom_Crop_Margin-1+y            ][2]=0;	
		if( Right__Crop_Margin > 0 )	for( y=0; y<Bimage->ysize+1; y+=2 ) Rimage->image[Left___Crop_Margin  +Bimage->xsize][Bottom_Crop_Margin-1+y            ][0]=Rimage->image[Left___Crop_Margin  +Bimage->xsize][Bottom_Crop_Margin-1+y            ][1]=Rimage->image[Left___Crop_Margin  +Bimage->xsize][Bottom_Crop_Margin-1+y            ][2]=0;		

		    return(TRUE );		
	 } else return(FALSE);
}


// Aligns the current image by shifting it sideways to match the previous image as much as possible. Required if there is significant drift from one image to the next. 

void align_images(int level) 
{
	char filename[STRING_LENGTH]; 
	int i, j, n, x, y; 
	int min_i, min_j, new_min_i, new_min_j, ai, aj;
	int search_i_min, search_i_max, search_j_min, search_j_max;
	float f1, f2; 
	float fractional_min_i, fractional_min_j; 
	float decimal_shift_x_px, decimal_shift_y_px;
	float difference, min_difference; 
	float d1, d3, e1, e3; 
	float difference_matrix[2*SHIFTXY+1][2*SHIFTXY+1];

	if( level == 0 ){				// Not shifting the first image, but just copying the Previous_Expanded_Image[][] over to Previous_Expanded_Image[][] in anticipation of the next level
		for( y=0; y<MAXY; y++ )
			for( x=0; x<MAXX; x++ )
				Previous_Expanded_Image[x][y] = Expanded_Image[x][y]; 
		Image_Shift_um_X[0] = Image_Shift_um_Y[0] = 0; 
		printf("  0.00      0.00    ");
		fprintf(CMDWindowFile, "0.000,0.000,");
		
	} else {

		for( j=-SHIFTXY; j<=SHIFTXY; j++ )
			for( i=-SHIFTXY; i<=SHIFTXY; i++ )
				difference_matrix[SHIFTXY+i][SHIFTXY+j] = BIG;				// Filling search matrix with large, positive values 

		min_i = new_min_i = search_i_min = search_i_max = Previous_Image_Shift_X_px;						// Starting search at previous best fit to save time
		min_j = new_min_j = search_j_min = search_j_max = Previous_Image_Shift_Y_px;	
//		min_i = min_j = new_min_i = new_min_j = 0;															// Starting search at zero image shift 
		min_difference = difference_matrix[SHIFTXY+min_i][SHIFTXY+min_j] = average_abs_grey_level_difference(min_i, min_j);		// Finding grey level difference in starting point separately 

		do{																	// Checking if neighbouring shifts give smaller difference and iterating until current shift gives local minimum in the difference_matrix[][]. Avoiding repeating the time consuming calculation in average_abs_grey_level_difference(). 
			min_i = new_min_i;	
			min_j = new_min_j;
			// Searching eight neighbouring pixels with increasing separation
			for( n=1; n<=5; n+=2 ) {
				if( difference_matrix[SHIFTXY+min_i+n][SHIFTXY+min_j  ] > BIG-1 ) { difference = difference_matrix[SHIFTXY+min_i+n][SHIFTXY+min_j  ] = average_abs_grey_level_difference(min_i+n, min_j  ); if( min_i+n > search_i_max ) search_i_max = min_i+n;													  }		if( difference < min_difference ) { new_min_i = min_i+n; new_min_j = min_j  ; min_difference = difference; }
				if( difference_matrix[SHIFTXY+min_i+n][SHIFTXY+min_j+n] > BIG-1 ) { difference = difference_matrix[SHIFTXY+min_i+n][SHIFTXY+min_j+n] = average_abs_grey_level_difference(min_i+n, min_j+n); if( min_i+n > search_i_max ) search_i_max = min_i+n; if( min_j+n > search_j_max ) search_j_max = min_j+n; }		if( difference < min_difference ) { new_min_i = min_i+n; new_min_j = min_j+n; min_difference = difference; }
				if( difference_matrix[SHIFTXY+min_i  ][SHIFTXY+min_j+n] > BIG-1 ) { difference = difference_matrix[SHIFTXY+min_i  ][SHIFTXY+min_j+n] = average_abs_grey_level_difference(min_i  , min_j+n);														 if( min_j+n > search_j_max ) search_j_max = min_j+n; }		if( difference < min_difference ) { new_min_i = min_i  ; new_min_j = min_j+n; min_difference = difference; }
				if( difference_matrix[SHIFTXY+min_i-n][SHIFTXY+min_j+n] > BIG-1 ) { difference = difference_matrix[SHIFTXY+min_i-n][SHIFTXY+min_j+n] = average_abs_grey_level_difference(min_i-n, min_j+n); if( min_i-n < search_i_min ) search_i_min = min_i-n; if( min_j+n > search_j_max ) search_j_max = min_j+n; }		if( difference < min_difference ) { new_min_i = min_i-n; new_min_j = min_j+n; min_difference = difference; }
				if( difference_matrix[SHIFTXY+min_i-n][SHIFTXY+min_j  ] > BIG-1 ) { difference = difference_matrix[SHIFTXY+min_i-n][SHIFTXY+min_j  ] = average_abs_grey_level_difference(min_i-n, min_j  ); if( min_i-n < search_i_min ) search_i_min = min_i-n;													  }		if( difference < min_difference ) { new_min_i = min_i-n; new_min_j = min_j  ; min_difference = difference; }
				if( difference_matrix[SHIFTXY+min_i-n][SHIFTXY+min_j-n] > BIG-1 ) { difference = difference_matrix[SHIFTXY+min_i-n][SHIFTXY+min_j-n] = average_abs_grey_level_difference(min_i-n, min_j-n); if( min_i-n < search_i_min ) search_i_min = min_i-n; if( min_j-n < search_j_min ) search_j_min = min_j-n; }		if( difference < min_difference ) { new_min_i = min_i-n; new_min_j = min_j-n; min_difference = difference; }
				if( difference_matrix[SHIFTXY+min_i  ][SHIFTXY+min_j-n] > BIG-1 ) { difference = difference_matrix[SHIFTXY+min_i  ][SHIFTXY+min_j-n] = average_abs_grey_level_difference(min_i  , min_j-n);														 if( min_j-n < search_j_min ) search_j_min = min_j-n; }		if( difference < min_difference ) { new_min_i = min_i  ; new_min_j = min_j-n; min_difference = difference; }
				if( difference_matrix[SHIFTXY+min_i+n][SHIFTXY+min_j-n] > BIG-1 ) { difference = difference_matrix[SHIFTXY+min_i+n][SHIFTXY+min_j-n] = average_abs_grey_level_difference(min_i+n, min_j-n); if( min_i+n > search_i_max ) search_i_max = min_i+n; if( min_j-n < search_j_min ) search_j_min = min_j-n; }		if( difference < min_difference ) { new_min_i = min_i+n; new_min_j = min_j-n; min_difference = difference; }
			}
		} while( new_min_i != min_i || new_min_j != min_j );

		d1 = difference_matrix[SHIFTXY+min_i-1][SHIFTXY+min_j  ] - difference_matrix[SHIFTXY+min_i][SHIFTXY+min_j];
		d3 = difference_matrix[SHIFTXY+min_i+1][SHIFTXY+min_j  ] - difference_matrix[SHIFTXY+min_i][SHIFTXY+min_j];
		fractional_min_i = (float)(0.5*(d1-d3)/(d1+d3));																// Finding bottom point of parabola through three points centered on lowest point in x-direction 

		e1 = difference_matrix[SHIFTXY+min_i  ][SHIFTXY+min_j-1] - difference_matrix[SHIFTXY+min_i][SHIFTXY+min_j];	
		e3 = difference_matrix[SHIFTXY+min_i  ][SHIFTXY+min_j+1] - difference_matrix[SHIFTXY+min_i][SHIFTXY+min_j];
		fractional_min_j = (float)(0.5*(e1-e3)/(e1+e3));																// Finding bottom point of parabola through three points centered on lowest point in y-direction 

		Previous_Image_Shift_X_px = min_i; 
		Previous_Image_Shift_Y_px = min_j; 

		if( fractional_min_i < 0 ) { fractional_min_i += 1; min_i -= 1; }			// Wanting only positive fractional_min values in the interval [0, 1> for easy interpolation below 
		if( fractional_min_j < 0 ) { fractional_min_j += 1; min_j -= 1; }			

		decimal_shift_x_px = min_i + fractional_min_i;
		decimal_shift_y_px = min_j + fractional_min_j;

		Image_Shift_um_X[level] = decimal_shift_x_px*Pixel_Size;
		Image_Shift_um_Y[level] = decimal_shift_y_px*Pixel_Size;

 		printf("%6.2f    %6.2f    ", Image_Shift_um_X[level], Image_Shift_um_Y[level]); 
 		fprintf(CMDWindowFile, "%5.3f,%5.3f,", Image_Shift_um_X[level], Image_Shift_um_Y[level]); 
	
		// Writing difference_matrix[][] to file for checking 
		fprintf(ShiftFile, "Image number %d,,,,,Shift X = %5.3f pixels,,,,,Shift Y = %5.3f pixels\nY \\ X,", level, decimal_shift_x_px, decimal_shift_y_px);		 
		for( i=search_i_min; i<=search_i_max; i++ ) fprintf(ShiftFile, "%d,", i);	fprintf(ShiftFile, "\n");
		for( j=search_j_min; j<=search_j_max; j++ ){
			fprintf(ShiftFile, "%d,", j); 
			for( i=search_i_min; i<=search_i_max; i++ ){
				if( (difference = difference_matrix[SHIFTXY+i][SHIFTXY+j]) < BIG-1 ) fprintf(ShiftFile, "%7.5f,", difference);
				else fprintf(ShiftFile, ",");
			}
			fprintf(ShiftFile, "\n");		
		}
		fprintf(ShiftFile, "\n");		// Blank line before next matrix 

		ai = abs(min_i);				// Shorthand since used often below 
		aj = abs(min_j);

		// Remapping the new image with correct fractional shift to Previous_Expanded_Image[][] since need to store it there anyway for shifting the next level image and since temporary storage is needed
		// In high contrast regions with alternating bright and dark pixels, this remapping will reduce the contrast and act as a kind of filter 
		for( y=-aj; y<Bimage->ysize+aj; y++ )
			for( x=-ai; x<Bimage->xsize+ai; x++ ){
				if( fractional_min_i + fractional_min_i < 1 )	f1 = (Expanded_Image[SHIFTXY+x+min_i+1][SHIFTXY+y+min_j  ]-Expanded_Image[SHIFTXY+x+min_i  ][SHIFTXY+y+min_j  ])*(  fractional_min_i) + (Expanded_Image[SHIFTXY+x+min_i  ][SHIFTXY+y+min_j+1]-Expanded_Image[SHIFTXY+x+min_i  ][SHIFTXY+y+min_j  ])*(  fractional_min_j) + Expanded_Image[SHIFTXY+x+min_i  ][SHIFTXY+y+min_j  ];
				else											f1 = (Expanded_Image[SHIFTXY+x+min_i  ][SHIFTXY+y+min_j+1]-Expanded_Image[SHIFTXY+x+min_i+1][SHIFTXY+y+min_j+1])*(1-fractional_min_i) + (Expanded_Image[SHIFTXY+x+min_i+1][SHIFTXY+y+min_j  ]-Expanded_Image[SHIFTXY+x+min_i+1][SHIFTXY+y+min_j+1])*(1-fractional_min_j) + Expanded_Image[SHIFTXY+x+min_i+1][SHIFTXY+y+min_j+1];
				if( fractional_min_j < fractional_min_i     )	f2 = (Expanded_Image[SHIFTXY+x+min_i  ][SHIFTXY+y+min_j  ]-Expanded_Image[SHIFTXY+x+min_i+1][SHIFTXY+y+min_j  ])*(1-fractional_min_i) + (Expanded_Image[SHIFTXY+x+min_i+1][SHIFTXY+y+min_j+1]-Expanded_Image[SHIFTXY+x+min_i+1][SHIFTXY+y+min_j  ])*(  fractional_min_j) + Expanded_Image[SHIFTXY+x+min_i+1][SHIFTXY+y+min_j  ];
				else											f2 = (Expanded_Image[SHIFTXY+x+min_i+1][SHIFTXY+y+min_j+1]-Expanded_Image[SHIFTXY+x+min_i  ][SHIFTXY+y+min_j+1])*(  fractional_min_i) + (Expanded_Image[SHIFTXY+x+min_i  ][SHIFTXY+y+min_j  ]-Expanded_Image[SHIFTXY+x+min_i  ][SHIFTXY+y+min_j+1])*(1-fractional_min_j) + Expanded_Image[SHIFTXY+x+min_i  ][SHIFTXY+y+min_j+1];
				Previous_Expanded_Image[SHIFTXY+x][SHIFTXY+y] = (float)(0.5*( f1 + f2 )); 	
			} 

		// Copying back the correctly shifted Expanded_Image[][] for further processing 
			for( y=0; y<MAXY; y++ )
				for( x=0; x<MAXX; x++ )
					Expanded_Image[x][y] = Previous_Expanded_Image[x][y];

		// Writing the shifted image also into the original Rimage (within the crop margins) to ensure that color dots will still match the displayed image 
		for( y=0; y<Bimage->ysize; y++ )
			for( x=0; x<Bimage->xsize; x++ )
				Rimage->image[x+Left___Crop_Margin][y+Bottom_Crop_Margin][0]=Rimage->image[x+Left___Crop_Margin][y+Bottom_Crop_Margin][1]=Rimage->image[x+Left___Crop_Margin][y+Bottom_Crop_Margin][2]=(unsigned char)(Expanded_Image[SHIFTXY+x][SHIFTXY+y]);	   
	}	

	// Writing SHIFTED image to file for checking. If only writing valid track images is turned on, also writing shifted images here and later deleting them after reading back in and adding overlays for videos 
	if( Write_1_SH_img || Write_7_VT_img || Make__7_VT_vid) {
		copy_bitmap(Rimage, Oimage);													// Copying original image into output image 
		for( y=0; y<Bimage->ysize; y++ )
			for( x=0; x<Bimage->xsize; x++ ) 
				Oimage->image[x+Left___Crop_Margin][y+Bottom_Crop_Margin][0] = (unsigned char)(Expanded_Image[SHIFTXY+x][SHIFTXY+y]); 	
		sprintf(filename, "%s/%s_01_%04d_Shifted.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*level);
		out_bitmap(Oimage, filename);
	}
}


// Since all the shifted images are read back at the end of execution to draw the valid tracks, using this function to draw dummy, non-shifted images when image alignment is turned off. 

void just_write_dummy_shifted_images(int level) 
{
	char filename[STRING_LENGTH]; 
	int x, y; 

	// Writing NON-SHIFTED image to file since needed if Write_7_VT_img or Make__7_VT_vid true for reading in fresh images.
	if( Write_7_VT_img || Make__7_VT_vid) {
		copy_bitmap(Rimage, Oimage);													// Copying original image into output image 
		for( y=0; y<Bimage->ysize; y++ )
			for( x=0; x<Bimage->xsize; x++ ) 
				Oimage->image[x+Left___Crop_Margin][y+Bottom_Crop_Margin][0] = (unsigned char)(Expanded_Image[SHIFTXY+x][SHIFTXY+y]); 	
		sprintf(filename, "%s/%s_01_%04d_Shifted.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*level);
		out_bitmap(Oimage, filename);
	}
	printf("   --        --     ");						// Writing blank info in shift columns on console
	fprintf(CMDWindowFile, "0,0,");						// Writing the same to Log window summary file
}


// Auxiliary function used often with different, specific shifts in align_images() above 

float average_abs_grey_level_difference(int i, int j)
{
	int x, y; 
	float diff_sum; 

	diff_sum = 0;
	for( y=0; y<Bimage->ysize; y++ )
		for( x=0; x<Bimage->xsize; x++ )
			diff_sum += abs( Previous_Expanded_Image[SHIFTXY+x][SHIFTXY+y] - Expanded_Image[SHIFTXY+x+i][SHIFTXY+y+j] );		 
	diff_sum /= Bimage->xsize*Bimage->ysize;		// Calculating the average difference per pixel to get a small number

	return(diff_sum);
}


// Runs a smoothing filter over Expanded_Image[][], automatically determines grey level limit, finds connected dark regions, identifies individual connected regions, scraps too small or large regions

void smooth_images_and_identify_cells(int level) 
{
	char filename[STRING_LENGTH], current_first_image_name[STRING_LENGTH*2], message[STRING_LENGTH], line[STRING_LENGTH]; 
	unsigned char color[3]; 
	short s; 
	int a, c, i, j, x, y, hist_index, top_i, left__i, right_i, center_i, grey_mode, histogram_nof__bins, histogram_half_bins, counter; 
	int box_height, box_width;
	int box_midpoint_y; 
	long sn, smoothed_histogram_only_background[501], smoothed_histogram_full_image[501], top_value; 
	float f, histogram_zoom_factor;
	double ssum, ssum2, smean, ssd, user_segmentation_grey_limit, lower_segmentation_grey_limit, upper_segmentation_grey_limit, cc, dd; 
	double log10_normalized_to_1_histogram[501]; 

	// Parabola fitting parameters 											
	int n, halflength;									 
	double aaa, ccc, eee, p;					// Coefficients for least square fit parabola to 2p+1 points in profile 
	double z0, z1, z2;							// Sum values for right hand side in matrix equation                    
	double A, B, C;								// Coefficients in z = Axx + Bx + C decribing the fitting parabolas		 

	histogram_half_bins = 100;
	histogram_nof__bins = 2*histogram_half_bins; 
	histogram_zoom_factor = (float)(histogram_nof__bins/10.0);		// The denominator gives the range of grey levels to include in the histogram (centred on the most common grey value)

	// For first image in each video and for first tuning image: Masking all bright and dark regions assumed to be cells and expanding the regions to be sure to mask everything but the background 
	sprintf(current_first_image_name, "%s/%s_02_%04d.bmp", InputFolderName, BMPfilename, First_Image_Number);
	if( level==0 && Tuning_Image_Code<2 ) { 

		if( FALSE && !strcmp(current_first_image_name, PreviousFirstImageName) ) {		// TEMPORARILY NOT DOING THIS TEST SO ALL RUNS WILL FIND OPTIMAL GREYLEVEL IN CASE ONLY THE Segmentation_Limit HAS BEEN CHANGED MANUALLY IN THE SETTINGS FILE
			printf("\nThe image %s \nhas been analyzed earlier and gave an optimal grey level limit of %5.3f. No need to do it again.\n", PreviousFirstImageName, Previous_Segmentation_Grey_Limit_From_First_Image);
		} else {

			// Running a small filter of local pixel and two nearest neighbours levels from Expanded_Image[][] to Smoothed_Image[][] in order to get a good spread of floating point grey level values for the histograms below 
			for( y=0; y<Bimage->ysize; y++ )
				for( x=0; x<Bimage->xsize; x++ ) 
	//			Smoothed_Image[SHIFTXY+x][SHIFTXY+y] = Expanded_Image[SHIFTXY+x][SHIFTXY+y];
				Smoothed_Image[SHIFTXY+x][SHIFTXY+y] = 0.2*Expanded_Image[SHIFTXY+x][SHIFTXY+y] + 0.15*(Expanded_Image[SHIFTXY+x+1][SHIFTXY+y]+Expanded_Image[SHIFTXY+x][SHIFTXY+y+1]+Expanded_Image[SHIFTXY+x-1][SHIFTXY+y]+Expanded_Image[SHIFTXY+x][SHIFTXY+y-1]) +  0.05*(Expanded_Image[SHIFTXY+x+1][SHIFTXY+y+1]+Expanded_Image[SHIFTXY+x-1][SHIFTXY+y+1]+Expanded_Image[SHIFTXY+x-1][SHIFTXY+y-1]+Expanded_Image[SHIFTXY+x+1][SHIFTXY+y-1]);

			// Calculating the mean and standard deviation of the Smoothed_Image[][]. 
			ssum = ssum2 = sn = 0; 
			for( y=0; y<Bimage->ysize; y++ )
				for( x=0; x<Bimage->xsize; x++ ) 
					if( (f=Smoothed_Image[SHIFTXY+x][SHIFTXY+y]) > 1 && f < 254 ) {				// Ignoring saturated pixels in black or white so including only greylevels between 1 and 254
						ssum  += f;
						ssum2 += f*f;
						sn++; 
					}
			if( sn > 1 ) {
				smean = ssum/sn;												// Average grey level in image, given as float 
				ssd   = sqrt( (ssum2-sn*smean*smean)/(sn-1) );					// Standard deviation of grey level 
			} else  { sprintf(message, "The smoothed image contains only white and black pixels. Please check data. \nExiting to avoid division by zero error when calculating mean and SD of graylevel in smoothed image.\n");  write_error_message(ErrorMessageFile, ERROR, message); } 	

			lower_segmentation_grey_limit = smean - Segmentation_Limit*ssd;		// New 2022-112-01: Using Segmentation_Limit here (also) to set limits for initial discrimination of cells from background. Tried since cells in collagen had much less smooth background 
			upper_segmentation_grey_limit = smean + Segmentation_Limit*ssd;		// Earlier comment: Seems to work well with mean ± one standard deviation as limits for dark and bright image regions 
			for( y=0; y<Bimage->ysize; y++ )
				for( x=0; x<Bimage->xsize; x++ ) 
					if(      Smoothed_Image[SHIFTXY+x][SHIFTXY+y] < lower_segmentation_grey_limit ) Segmented_Image[SHIFTXY+x][SHIFTXY+y] =   1;		// NB: Coding segmented pixel as 1 to distinquish from NEGATIVE numbered cells below 
					else if( Smoothed_Image[SHIFTXY+x][SHIFTXY+y] > upper_segmentation_grey_limit ) Segmented_Image[SHIFTXY+x][SHIFTXY+y] =   1;		// Masking also pixels that are BRIGHTER than normal in order to eliminate those from the background 
					else																			Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 255;		// Temporarily coding background as 255 to distinquish holes in cells from the continuous background after filling it with 0 below 

			// Dilating the the smoothed+segmented image several times to grow mask beyond all cells and leaving only background outside mask for grey level analysis, using Tuning_Image[][] as temporary storage 
			for( i=0; i<2; i++ ){				
				for(x=0; x<Bimage->xsize; x++) 		
					for(y=0; y<Bimage->ysize; y++)
						if( Segmented_Image[SHIFTXY+x][SHIFTXY+y] == 255 ){															// Pixel is in background 
								 if( Segmented_Image[SHIFTXY+x+1][SHIFTXY+y  ] == 1 ) Tuning_Image[SHIFTXY+x][SHIFTXY+y] = 1; 
							else if( Segmented_Image[SHIFTXY+x-1][SHIFTXY+y  ] == 1 ) Tuning_Image[SHIFTXY+x][SHIFTXY+y] = 1; 
							else if( Segmented_Image[SHIFTXY+x  ][SHIFTXY+y+1] == 1 ) Tuning_Image[SHIFTXY+x][SHIFTXY+y] = 1;		// Setting to masked region if four-connected to mask 
							else if( Segmented_Image[SHIFTXY+x  ][SHIFTXY+y-1] == 1 ) Tuning_Image[SHIFTXY+x][SHIFTXY+y] = 1; 
		/*					else if( Segmented_Image[SHIFTXY+x+1][SHIFTXY+y+1] == 1 ) Tuning_Image[SHIFTXY+x][SHIFTXY+y] = 1; 
							else if( Segmented_Image[SHIFTXY+x-1][SHIFTXY+y-1] == 1 ) Tuning_Image[SHIFTXY+x][SHIFTXY+y] = 1; 
							else if( Segmented_Image[SHIFTXY+x-1][SHIFTXY+y+1] == 1 ) Tuning_Image[SHIFTXY+x][SHIFTXY+y] = 1;		 
							else if( Segmented_Image[SHIFTXY+x+1][SHIFTXY+y-1] == 1 ) Tuning_Image[SHIFTXY+x][SHIFTXY+y] = 1;		*/ 
							else													  Tuning_Image[SHIFTXY+x][SHIFTXY+y] = 255;
						} else														  Tuning_Image[SHIFTXY+x][SHIFTXY+y] = 1;		// Copying over masked pixels
				for(x=0; x<Bimage->xsize; x++) 		
					for(y=0; y<Bimage->ysize; y++)
						if( Tuning_Image[SHIFTXY+x][SHIFTXY+y] == 255 ){															// Pixel is in background 
								 if( Tuning_Image[SHIFTXY+x+1][SHIFTXY+y  ] == 1 ) Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 1; 
							else if( Tuning_Image[SHIFTXY+x-1][SHIFTXY+y  ] == 1 ) Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 1; 
							else if( Tuning_Image[SHIFTXY+x  ][SHIFTXY+y+1] == 1 ) Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 1;		// Setting to masked region if four-connected to mask  
							else if( Tuning_Image[SHIFTXY+x  ][SHIFTXY+y-1] == 1 ) Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 1; 
		/*					else if( Tuning_Image[SHIFTXY+x+1][SHIFTXY+y+1] == 1 ) Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 1;		// 8 connected option made regions very rectangular 
							else if( Tuning_Image[SHIFTXY+x-1][SHIFTXY+y-1] == 1 ) Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 1;	
							else if( Tuning_Image[SHIFTXY+x-1][SHIFTXY+y+1] == 1 ) Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 1;		 
							else if( Tuning_Image[SHIFTXY+x+1][SHIFTXY+y-1] == 1 ) Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 1; 
							else													   Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 255;		*/ 
						} else														   Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 1;		// Copying back masked pixels
			}

			// Writing FIRST IMAGE WITH CELLS COVERED IN BLACK to file for checking		
			copy_bitmap(Rimage, Oimage);	 
			for( y=0; y<Bimage->ysize; y++ )
				for( x=0; x<Bimage->xsize; x++ ) 
					if( Segmented_Image[SHIFTXY+x][SHIFTXY+y] == 1 ) Oimage->image[x+Left___Crop_Margin][y+Bottom_Crop_Margin][0] = Oimage->image[x+Left___Crop_Margin][y+Bottom_Crop_Margin][1] = Oimage->image[x+Left___Crop_Margin][y+Bottom_Crop_Margin][2] = (unsigned char)0;			// (unsigned char)(Smoothed_Image[SHIFTXY+x][SHIFTXY+y]); 
					else  											 Oimage->image[x+Left___Crop_Margin][y+Bottom_Crop_Margin][0] = Oimage->image[x+Left___Crop_Margin][y+Bottom_Crop_Margin][1] = Oimage->image[x+Left___Crop_Margin][y+Bottom_Crop_Margin][2] = Bimage->image[x][y][0];	// Showing original image in non-masked regions for checking 
			sprintf(filename, "%s/%s____%04d_First_Image_Background_%02d_iterations.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*level, i);
			out_bitmap(Oimage, filename);		//*/
		} 
	} // End of FIRST "	if( level==0 && Tuning_Image_Code<2 )"

	// Gaussian smooting of Expanded_Image[][] which copies it over to the Smoothed_Image[][] array used below for more precise segmentation
	gaussian_smooth(Gaussian_Filter_Halfwidth);			// CONSIDER LATER TO LET FILTER ROUTINE ALSO CALCULATE MEAN AND SD AND RETURN THESE PARAMETERS

	// Writing SMOOTHED image to file for checking 
	if( Write_2_SM_img ) {
		copy_bitmap(Rimage, Oimage);													// Copying original image into output image 
		for( y=0; y<Bimage->ysize; y++ )
			for( x=0; x<Bimage->xsize; x++ ) 
				Oimage->image[x+Left___Crop_Margin][y+Bottom_Crop_Margin][0] = Oimage->image[x+Left___Crop_Margin][y+Bottom_Crop_Margin][1] = Oimage->image[x+Left___Crop_Margin][y+Bottom_Crop_Margin][2] = (unsigned char)(Smoothed_Image[SHIFTXY+x][SHIFTXY+y]); 	
		if( Add_Scalebar ) draw_scalebar(); 
		sprintf(filename, "%s/%s_02_%04d_Smoothed.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*level);
		out_bitmap(Oimage, filename);
	}

	// Calculating the mean and standard deviation of the Smoothed_Image[][]
	ssum = ssum2 = sn = 0; 
	for( y=0; y<Bimage->ysize; y++ )
		for( x=0; x<Bimage->xsize; x++ ) 
			if( (f=Smoothed_Image[SHIFTXY+x][SHIFTXY+y]) > 1 && f < 254 ) {		// Ignoring saturated pixels in black or white so including only greylevels between 1 and 254
				ssum  += f;
				ssum2 += f*f;
				sn++; 
			}
	if( sn > 1 ) {
		smean = ssum/sn;														// Average grey level in image, given as float 
		ssd   = sqrt( (ssum2-sn*smean*smean)/(sn-1) );							// Standard deviation of grey level 
		user_segmentation_grey_limit = smean-Segmentation_Limit*ssd;			// Calculating and displaying SD based grey level to compare with automatically found level based on first image below 	
	} else  { sprintf(message, "The smoothed image contains only white and black pixels. Please check data. \nExiting to avoid division by zero error when calculating mean and SD of graylevel in smoothed image.\n");  write_error_message(ErrorMessageFile, ERROR, message); } 	

	// For first image in each video and for first tuning image: Calculating a detailed grey level histogram of the smoothed background (which will roughly have a gaussian shape), fitting a parabola to the logarithm of the normalized histogram and then finding the best grey level limit for segmentation 
	if( level==0 && Tuning_Image_Code<2 ) {			

		if( FALSE && !strcmp(current_first_image_name, PreviousFirstImageName) ) { 		// TEMPORARILY NOT DOING THIS TEST SO ALL RUNS WILL FIND OPTIMAL GREYLEVEL IN CASE THE Segmentation_Limit HAS BEEN CHANGED MANUALLY IN THE SETTINGS FILE
			Segmentation_Grey_Limit_From_First_Image = Previous_Segmentation_Grey_Limit_From_First_Image; 
		} else {
			strcpy(PreviousFirstImageName, current_first_image_name);		// Making sure to update the PreviousFirstImageName when a new video is being processed to prevent the wrong segmentation limit to be used in the next run 

			// FIRST finding the most common grey value in the smoothed image (=grey_mode) by generating a histogram over all grey values that are not masked, including all grey values from 0 to 256
			for( i=0; i<256; i++ ) smoothed_histogram_only_background[i] = 0;			// Just making sure that the histogram is all zero before starting 
			for( y=0; y<Bimage->ysize; y++ )
				for( x=0; x<Bimage->xsize; x++ ) {
						hist_index = (int)( Smoothed_Image[SHIFTXY+x][SHIFTXY+y] + 0.5 );		
						if( hist_index <   0 ) hist_index =   0;					
						if( hist_index > 255 ) hist_index = 255;					
						if( Segmented_Image[SHIFTXY+x][SHIFTXY+y] == 255 ) smoothed_histogram_only_background[hist_index]++;					// The pixel is not masked as near a cell, so included in the background histogram 
					}
			top_value = grey_mode = 0; 
			for( i=1; i<255; i++ )														// Finding highest value while excluding edge points to make sure top_i is never there
				if( smoothed_histogram_only_background[i] > top_value ) { top_value = smoothed_histogram_only_background[i]; grey_mode = i; }
			// Now grey_mode is a robust estimate of the most common background grey level to the nearest integer  

			// THEN finding a more precice value of the most common grey value in the smoothed image by generating a histogram over a narrower range of grey values that are not masked, and fitting a parabola around the peak of the log10 of the normalized histogram
			for( i=0; i<histogram_nof__bins; i++ ) smoothed_histogram_full_image[i] = smoothed_histogram_only_background[i] = 0;				// Clearing the histogram arrays  
			for( y=0; y<Bimage->ysize; y++ )
				for( x=0; x<Bimage->xsize; x++ ) {
						hist_index = (int)( histogram_zoom_factor*(Smoothed_Image[SHIFTXY+x][SHIFTXY+y]-grey_mode) + histogram_half_bins );		// Zooming in the histogram around the most common grey value (grey_mode) found in the first iteration. 
						if( hist_index <                     0 ) hist_index = 0;					
						if( hist_index > histogram_nof__bins-1 ) hist_index = histogram_nof__bins-1;					
						smoothed_histogram_full_image[hist_index]++; 
						if( Segmented_Image[SHIFTXY+x][SHIFTXY+y] == 255 ) smoothed_histogram_only_background[hist_index]++;					// The pixel is not masked as near a cell, so included in the background histogram 
					}
			top_value = top_i = 0; 
			for( i=1; i<histogram_nof__bins-1; i++ )														// Finding highest value while excluding edge points to make sure top_i is never there
				if( smoothed_histogram_only_background[i] > top_value ) { top_value = smoothed_histogram_only_background[i]; top_i = i; }
			for( i=0; i<histogram_nof__bins; i++ )															// Normalizing and taking the logarithm to convert a gaussian to a parabola
				if( smoothed_histogram_only_background[i] > 0 )	log10_normalized_to_1_histogram[i] = log10(1.0*smoothed_histogram_only_background[i]/top_value);		// Normalizing so max value is 1 before taking lg(), after which the top point is 0, the 10 % level is -1, the 1 % level is -2 and the 0.1 % level is -3, etc. 
				else									        log10_normalized_to_1_histogram[i] = -10;																// Flagging zero counts as a large negative value which should not occur otherwise 

			left__i = top_i;	while( log10_normalized_to_1_histogram[--left__i] > -1.5 );		// Taking -1.5 as left  limit, corresponding to 1/30th of the max value
			right_i = top_i;	while( log10_normalized_to_1_histogram[++right_i] > -0.5 ); 	// Taking -0.5 as right limit, corresponding to 1/3rd  of the max value since stronger non-gaussian tail at right

			// Fitting a parabola to the points above the log10_limit value  
			halflength = (int)( 0.5*(right_i-left__i) );			// Handling asymmetric limit points 
			center_i   = (int)( 0.5*(right_i+left__i) );			// Handling asymmetric limit points 
			p = (double)halflength;									// Shorthand since used a lot below 
			aaa = 2*(p*p*p*p*p/5 + p*p*p*p/2 + p*p*p/3 - p/30);		// Calculating the matrix coefficients 
			ccc = p*(p+1)*(2*p+1)/3;
			eee = 2*p+1;
			z0=0; for( n=-halflength; n<=halflength; n++ ) z0 +=     log10_normalized_to_1_histogram[center_i+n];
			z1=0; for( n=-halflength; n<=halflength; n++ ) z1 +=   n*log10_normalized_to_1_histogram[center_i+n]; 
			z2=0; for( n=-halflength; n<=halflength; n++ ) z2 += n*n*log10_normalized_to_1_histogram[center_i+n]; 
			A = ( ccc*z0 - eee*z2 )/( ccc*ccc - aaa*eee );	
			B = z1/ccc;												// Coefficients in fitted parabola 
			C = ( ccc*z2 - aaa*z0 )/( ccc*ccc - aaa*eee );		

	//		cc = B*B - 4.0*A*(C+2.5);	// +2.5 corresponds to a count at this greylevel which is 0.3 % of the max value in the nearly gaussian grey level distribution in the filtered background (+1 would be 10 %, +2 would be 1 %, +3 would be 0.1 %, etc.) 
			cc = B*B - 4.0*A*(C+3.0);	//   +3 corresponds to a count at this greylevel which is 0.1 % of the max value in the nearly gaussian grey level distribution in the filtered background , thus selecting a grey level limit where only 1/1000 of the background pixels are darker
			if( cc>0 && A<0 ) { 
				dd = 0.5*(-B + sqrt(cc))/A;		
				Segmentation_Grey_Limit_From_First_Image = (float)(grey_mode+(center_i+dd-histogram_half_bins)/histogram_zoom_factor);					
			} else { 
				printf("\nWarning: Could not fit parabola to top of smoothed grey image histogram.\nInstead using the SD based segmentation limit (%4.2lf)\n", user_segmentation_grey_limit); 
				Segmentation_Grey_Limit_From_First_Image = (float)user_segmentation_grey_limit;
			} 
			if( Use_SD_Based_Segmentation_Grey_Limit ) Segmentation_Grey_Limit_From_First_Image = (float)user_segmentation_grey_limit;			// Option to override the auto segmentation limit if the parameter Segmentation limit is given as a negative number in celltraxx_settings.txt

			for( i=0; i<histogram_nof__bins; i++ )					// Writing histograms and fitted parabola to file for possible checking
				if( log10_normalized_to_1_histogram[i] == -10 )		fprintf(GreyHistogramsFile,	"%s%04d.bmp,%d,%5.3lf,%ld,%ld,,%5.3lf\n",       BMPfilename, First_Image_Number, i, grey_mode+(i-histogram_half_bins)/histogram_zoom_factor, smoothed_histogram_full_image[i], smoothed_histogram_only_background[i],                                     A*(i-center_i)*(i-center_i) + B*(i-center_i) + C);	// Skipping invalid log values for zero count grey levels
				else												fprintf(GreyHistogramsFile,	"%s%04d.bmp,%d,%5.3lf,%ld,%ld,%5.3lf,%5.3lf\n", BMPfilename, First_Image_Number, i, grey_mode+(i-histogram_half_bins)/histogram_zoom_factor, smoothed_histogram_full_image[i], smoothed_histogram_only_background[i], log10_normalized_to_1_histogram[i], A*(i-center_i)*(i-center_i) + B*(i-center_i) + C);	
			fprintf(GreyHistogramsFile, "\n");						// One blank line before next first image 

			// Writing updated info back to the tuning image file as a way to remember the segmentation limit if the same image is analysed again and if tuning the middle and last image of the same video.
			sprintf(filename, "C:/celltraxx_system/celltraxx_tuning_info.txt"); 
			if( (TuneFile=fopen(filename, "w"))==NULL ){ sprintf(message, "Can't open file celltraxx_tuning_info.txt for writing.\nPlease check if the file is open in another application.\n");  write_error_message(ErrorMessageFile, ERROR, message); exit(0); }
			fprintf(TuneFile, "First image grey level limit %5.3f\n", Segmentation_Grey_Limit_From_First_Image);
			fprintf(TuneFile, "First image path and name    %s\n",    PreviousFirstImageName);
			fclose(TuneFile);	//*/
		}
	} // End of SECOND "	if( level==0 && Tuning_Image_Code<2 )"

	if( Tuning_Image_Code > 1 ) { 		// In order to use the Segmentation_Grey_Limit_From_First_Image also when analyzing middle and last tuning image, reading the segmentation limit from disk instead of doing the above analyses. Since tuning always starts with the first selected image, the correct segmentation limit exists in the file. 
		sprintf(filename, "C:/celltraxx_system/celltraxx_tuning_info.txt"); 
		if( (TuneFile=fopen(filename, "r"))==NULL ){ sprintf(message, "Can't open file celltraxx_tuning_info.txt for reading.\nPlease check if the file is open in another application.\n");  write_error_message(ErrorMessageFile, ERROR, message); }
		if( (fgets(line, STRING_LENGTH, TuneFile)==NULL) || sscanf(line, "First image grey level limit %f", &Segmentation_Grey_Limit_From_First_Image)<1  || Segmentation_Grey_Limit_From_First_Image <0 ) { sprintf(message, "Missing or invalid value for 'First image grey level limit' (%5.3f) in file %s. Should be a positive number.", Segmentation_Grey_Limit_From_First_Image, filename); write_error_message(ErrorMessageFile, ERROR, message); }
		fclose(TuneFile);
	}

	if( Segmentation_Grey_Limit_From_First_Image > smean ) {			// If auto segmentation clearly fails, using instead the SD based value of greylevel_mean - Segmentation_Limit*greylevel_SD 
		printf("Warning: The automatic segmentation limit (%4.2f) is brighter than the mean grey level (%4.2lf).\nInstead using the SD based segmentation limit (%4.2lf)\n", Segmentation_Grey_Limit_From_First_Image, smean, user_segmentation_grey_limit); 
		Segmentation_Grey_Limit_From_First_Image = (float)user_segmentation_grey_limit;
	} 
	if( Use_SD_Based_Segmentation_Grey_Limit ) Segmentation_Grey_Limit_From_First_Image = (float)user_segmentation_grey_limit;

	printf("%6.2f  %5.2f      %6.2f         %6.2f      ", smean, ssd, user_segmentation_grey_limit, Segmentation_Grey_Limit_From_First_Image); 
	fprintf(CMDWindowFile, "%5.3f,%5.3f,%5.3f,%5.3f,",    smean, ssd, user_segmentation_grey_limit, Segmentation_Grey_Limit_From_First_Image); 

	// Making copy of the previous segmented image for use below to generate overlap image between current and previous to separate merged cells 
	if( level == 0 ) {
		for(x=0; x<Bimage->xsize; x++) 		
			for(y=0; y<Bimage->ysize; y++)
				Previous_Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 1;		// For zeroth level using no mask so all cells are accepted. At other levels, using Previous_Segmented_Image[][] copied from Segmented_Image[][] after scrapping too small and large regions near the end of this function 
	} /*else {
		for(x=0; x<Bimage->xsize; x++) 		
			for(y=0; y<Bimage->ysize; y++)
				Previous_Segmented_Image[SHIFTXY+x][SHIFTXY+y] = Segmented_Image[SHIFTXY+x][SHIFTXY+y];		
	} //*/
	
	// Segmenting pixels in the Smoothed_Image[][] that are darker than Segmentation_Grey_Limit_From_First_Image 
	for( y=0; y<Bimage->ysize; y++ )
		for( x=0; x<Bimage->xsize; x++ ) 
			if(   Smoothed_Image[SHIFTXY+x][SHIFTXY+y] < Segmentation_Grey_Limit_From_First_Image ) Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 1;		// NB: Coding segmented pixel as 1 to distinquish from NEGATIVE numbered cells below 
			else Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 255;														// Temporarily coding background as 255 to distinquish holes in cells from the continuous background after filling it with 0 below 

	// Writing SEGMENTED image to file for checking. Drawing only cell contours for more easy comparison.
	if( Write_3_SE_img ) {
		copy_bitmap(Rimage, Oimage);								// Copying original image into output image 
		convert2color(Oimage);										// To allow showing different cells with different colors 
		for( y=0; y<Bimage->ysize; y++ )
			for( x=0; x<Bimage->xsize; x++ ) 
				if( Segmented_Image[SHIFTXY+x][SHIFTXY+y]==1 && (Segmented_Image[SHIFTXY+x+1][SHIFTXY+y]!=1 || Segmented_Image[SHIFTXY+x-1][SHIFTXY+y]!=1 || Segmented_Image[SHIFTXY+x][SHIFTXY+y+1]!=1 || Segmented_Image[SHIFTXY+x][SHIFTXY+y-1]!=1 ) ) {
					cellnumber2monocolor(1,color,6);	// The colorcode sets channel by bit number 1=blue 2=green 4=red. Can be combined to give shades of cyan=3, magenta=5, yellow=6. 
					for( i=0; i<3; i++) Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][i] = color[i]; 			
				}
		if( Add_Scalebar ) draw_scalebar(); 
		sprintf(filename, "%s/%s_03_%04d_Segmented.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*level);
		out_bitmap(Oimage, filename);
	}
	
	Min_Diameter_px = (int)(Min_Diameter_um/Pixel_Size+0.5);					// Shorthand since used several places below 
	Max_Diameter_px = (int)(Max_Diameter_um/Pixel_Size+0.5);		 
	Mid_Diameter_px = (int)(0.5*(Min_Diameter_px+Max_Diameter_px));
	Min_Cell_Pixels = (int)(Pi*Min_Diameter_px*Min_Diameter_px/4);				// Limit for scrapping segmented features smaller than the circular area given by the smallest cell diameter.  
	Max_Cell_Pixels = (int)(Pi*Max_Diameter_px*Max_Diameter_px/4);				// Limit for scrapping segmented features larger  than the circular area given by the largest  cell diameter.
	Mid_Cell_Pixels = (int)(Pi*Mid_Diameter_px*Mid_Diameter_px/4);				// Limit for cutting cells at concave points.  

	// Painting background with zeros before filling enclosed holes (still marked by 255) inside cells. Filling from all four corners in case a dense row of cells near a corner would make a barrier against filling all background pixels
	x=i=0;												 
	while( i<3*Min_Diameter_px && x<MAXX )														// Lower left corner
		if( Segmented_Image[SHIFTXY+x++][SHIFTXY] == 255 ) i++; else i=0;						// Moving along bottom pixel row until having found a connected strip outside a cell with length three typical cell diameters. This should ensure that we are outside any hole in a cell. 
	j = fill_segmented(255, 0, x, 0);															// Filling all of background with 0, but leaving holes inside cells at 255 
	x=i=0; 	 
	while( i<3*Min_Diameter_px && x<MAXX )														// Upper left corner
		if( Segmented_Image[SHIFTXY+x++][SHIFTXY+Bimage->ysize-1] == 255 ) i++; else i=0;		// Moving along top pixel row until having found a connected strip outside a cell with length three typical cell diameters. This should ensure that we are outside any hole in a cell. 
	j = fill_segmented(255, 0, x, Bimage->ysize-1);												// Filling all of background with 0, but leaving holes inside cells at 255 
	x=Bimage->xsize-1; i=0;	 
	while( i<3*Min_Diameter_px && x>0 )															// Lower right corner
		if( Segmented_Image[SHIFTXY+x--][SHIFTXY] == 255 ) i++; else i=0;						// Moving along bottom pixel row until having found a connected strip outside a cell with length three typical cell diameters. This should ensure that we are outside any hole in a cell. 
	j = fill_segmented(255, 0, x, 0);															// Filling all of background with 0, but leaving holes inside cells at 255 
	x=Bimage->xsize-1; i=0;  
	while( i<3*Min_Diameter_px && x>0 )															// Upper right corner
		if( Segmented_Image[SHIFTXY+x--][SHIFTXY+Bimage->ysize-1] == 255 ) i++; else i=0;		// Moving along top pixel row until having found a connected strip outside a cell with length three typical cell diameters. This should ensure that we are outside any hole in a cell. 
	j = fill_segmented(255, 0, x, Bimage->ysize-1);												// Filling all of background with 0, but leaving holes inside cells at 255 

	// Filling holes: Going through the Segmented_Image[][] image and replacing any pixels of value 255 with 1 since these are internal holes in cells that were not affected when the background was filled in 
	for( y=0; y<Bimage->ysize; y++ )
		for( x=0; x<Bimage->xsize; x++ ) 
			if(  Segmented_Image[SHIFTXY+x][SHIFTXY+y] == 255 ) Segmented_Image[SHIFTXY+x][SHIFTXY+y] = 1;	// Now the Segmented_Image[][] has background = 0 and cells = 1 

	// Analysing large regions and cutting them at narrow and concave points, assuming this means that two cells have merged 
	Nof__Cells[level]=3;				// Here the temporary numbering starts at 4 in order to reserve lower number for overlap coding below 
	for(x=0; x<Bimage->xsize; x++) 		
		for(y=0; y<Bimage->ysize; y++)
			if( Segmented_Image[SHIFTXY+x][SHIFTXY+y]==1 ){									// Pixel belongs to region marked with NEGATIVE number above  
				i = ++Nof__Cells[level];
				if( i == MAXCELLS-1 ) { sprintf(message, "The number of cells (Nof__Cells[level] = %d) has reached the \nsize of the arrays which store cell data (dimension = %d) at image number %d.\nReduce the number of images or increase MAXCELLS and recompile celltraxx.c.\n", Nof__Cells[level], MAXCELLS, level+First_Image_Number);  write_error_message(ErrorMessageFile, ERROR, message); }
				a = fill_segmented_box(1, Nof__Cells[level], x, y, &Cells[level][Nof__Cells[level]].min_x, &Cells[level][Nof__Cells[level]].max_x, &Cells[level][Nof__Cells[level]].min_y, &Cells[level][Nof__Cells[level]].max_y);			// Filling cell with TEMPORARY POSITIVE identification number. Need to use fill_segmented_box() since cut_concave_region() reads the box coordinates Cells[level][Nof__Cells[level]].min_x, etc. 
				if( Cutting_Diameter_px>0 && a>Min_Cell_Pixels ) cut_concave_region(x, y, level);				// Special routine for cutting at convex, narrow passages. Cut pixels are coded as CUT in Segmented_Image[][]
			} 

	// Writing filled and CUT image to file for checking. Drawing only cell contours for more easy comparison.
	if( Write_4_CU_img ) {
		copy_bitmap(Rimage, Oimage);										// Copying original image into output image 
		convert2color(Oimage);												// To allow showing different cells with different colors 
		for( y=0; y<Bimage->ysize; y++ )
			for( x=0; x<Bimage->xsize; x++ )							//																				 Blue																			   Green																			   Red
				if(      Segmented_Image[SHIFTXY+x][SHIFTXY+y] == CUT  ) { Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0] = (unsigned char)255; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1] = (unsigned char)  0; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2] = (unsigned char)  0; }  	// Blue cutting lines		
				else if( Segmented_Image[SHIFTXY+x][SHIFTXY+y] == CUT2 ) { Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0] = (unsigned char)  0; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1] = (unsigned char)  0; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2] = (unsigned char)  0; } 	// Black all potential cut lines		
				else if( Segmented_Image[SHIFTXY+x][SHIFTXY+y] == FIT  ) { Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0] = (unsigned char)255; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1] = (unsigned char)  0; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2] = (unsigned char)255; } 	// Magenta for narrow concave points	
				else if( Segmented_Image[SHIFTXY+x][SHIFTXY+y] == FIT2 ) { Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0] = (unsigned char)255; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1] = (unsigned char)255; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2] = (unsigned char)  0; } 	// Cyan    for narrow convex points 		
				else if( Segmented_Image[SHIFTXY+x][SHIFTXY+y] == FIT3 ) { Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0] = (unsigned char)255; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1] = (unsigned char)255; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2] = (unsigned char)  0; } 	// Yellow  for narrow  flat  points		
				else if( (s=Segmented_Image[SHIFTXY+x][SHIFTXY+y])>0 && (Segmented_Image[SHIFTXY+x+1][SHIFTXY+y]!=s || Segmented_Image[SHIFTXY+x-1][SHIFTXY+y]!=s || Segmented_Image[SHIFTXY+x][SHIFTXY+y+1]!=s || Segmented_Image[SHIFTXY+x][SHIFTXY+y-1]!=s ) ) {
					                                                       Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0] = (unsigned char)  0; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1] = (unsigned char)255; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2] = (unsigned char)  0; 	// Green cell outlines		
				}
		if( Add_Scalebar ) draw_scalebar(); 
		sprintf(filename, "%s/%s_04_%04d_Cut.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*level);
		out_bitmap(Oimage, filename);	
	}
	// Now Segmented_Image[][] has background 0, cells filled with positive numbers 4, 5, 6,... and cut lines coded with CUT

	// Checking if pixels in Segmented_Image[][] overlap with the Previous_Segmented_Image[][] as a trick to separate regions which were apart in the previous image but are merged in the current image 
	for(x=0; x<Bimage->xsize; x++) 		
		for(y=0; y<Bimage->ysize; y++)
			if( Segmented_Image[SHIFTXY+x][SHIFTXY+y]>3 ) {
				if( Keep_Previous_Cells && Previous_Segmented_Image[SHIFTXY+x][SHIFTXY+y]!=0 ) Segmented_Image[SHIFTXY+x][SHIFTXY+y]=2;		// Coding as 2 if overlap with previous segmented image and this mode is active 
				else																		   Segmented_Image[SHIFTXY+x][SHIFTXY+y]=1;		// Coding as 1 if segmented in current image but no overlap with previous image
			} else if( Segmented_Image[SHIFTXY+x][SHIFTXY+y] == CUT ||  Segmented_Image[SHIFTXY+x][SHIFTXY+y] == FIT || Segmented_Image[SHIFTXY+x][SHIFTXY+y] == FIT2 ) Segmented_Image[SHIFTXY+x][SHIFTXY+y]=0;		// Erasing cutting lines since displayed in above image and not needed any more
	// Now Segmented_Image[][] has cells overlapping with previous image = 2, new cell regions = 1 and background = 0 

/*	if( Write_4_CU_img ) {
		copy_bitmap(Rimage, Oimage);								
		convert2color(Oimage);							
		for( y=0; y<Bimage->ysize; y++ )
			for( x=0; x<Bimage->xsize; x++ ) 
				     if( Segmented_Image[SHIFTXY+x][SHIFTXY+y]==2 ) { Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0]=Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1]=(unsigned char)0; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2]=(unsigned char)255; } 		
				else if( Segmented_Image[SHIFTXY+x][SHIFTXY+y]==1 ) { Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0]=(unsigned char)0; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1]=Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2]=(unsigned char)255; } 		
		if( Add_Scalebar ) draw_scalebar(); 
		sprintf(filename, "%s/%s_04_1_%04d_Overlap-regions.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*level);
		out_bitmap(Oimage, filename);	
	} //*/
	
	// Numbering OVERLAPPING segmented areas with NEGATIVE numbers to distinguish them below when giving them positive numbers. 
	counter = 1;																							
	if( Keep_Previous_Cells ) { 
		for(x=0; x<Bimage->xsize; x++) 		
			for(y=0; y<Bimage->ysize; y++)
				if( (s=Segmented_Image[SHIFTXY+x][SHIFTXY+y]) == 2 ){						// Pixel identified as dark region, potential part of cell which also overlaps previous image 
					a = fill_segmented(s, -counter, x, y);									// Filling cell with NEGATIVE number
					counter++;
				}

/*		if(  Write_4_CU_img ) {
			copy_bitmap(Rimage, Oimage);					
			convert2color(Oimage);								
			for( y=0; y<Bimage->ysize; y++ )
				for( x=0; x<Bimage->xsize; x++ ) 
					if( (s=Segmented_Image[SHIFTXY+x][SHIFTXY+y])!=0 ) {
						cellnumber2color(s,color);
						for( i=0; i<3; i++) Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][i] = color[i]; 			
					}
			if( Add_Scalebar ) draw_scalebar(); 
			sprintf(filename, "%s/%s_04_2_%04d_Overlapping-regions-filled.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*level);
			out_bitmap(Oimage, filename);	
		} //*/

		// Expanding overlap regions that were filled with identifying negative numbers until borders in current Segmented_Image[][] are reached. Using Tuning_Image[][] as intermediate storage.
		do{
			n=0; 
			for(x=0; x<Bimage->xsize; x++) 		
				for(y=0; y<Bimage->ysize; y++)
					if( Segmented_Image[SHIFTXY+x][SHIFTXY+y]==1 ) {
							 if( (s=Segmented_Image[SHIFTXY+x+1][SHIFTXY+y  ])<0 ) { Tuning_Image[SHIFTXY+x][SHIFTXY+y]=s; n++; }
						else if( (s=Segmented_Image[SHIFTXY+x-1][SHIFTXY+y  ])<0 ) { Tuning_Image[SHIFTXY+x][SHIFTXY+y]=s; n++; }
						else if( (s=Segmented_Image[SHIFTXY+x  ][SHIFTXY+y+1])<0 ) { Tuning_Image[SHIFTXY+x][SHIFTXY+y]=s; n++; }
						else if( (s=Segmented_Image[SHIFTXY+x  ][SHIFTXY+y-1])<0 ) { Tuning_Image[SHIFTXY+x][SHIFTXY+y]=s; n++; }
						else Tuning_Image[SHIFTXY+x][SHIFTXY+y] = Segmented_Image[SHIFTXY+x][SHIFTXY+y];
					} else   Tuning_Image[SHIFTXY+x][SHIFTXY+y] = Segmented_Image[SHIFTXY+x][SHIFTXY+y];
			for(x=0; x<Bimage->xsize; x++) 		
				for(y=0; y<Bimage->ysize; y++)
					if( Tuning_Image[SHIFTXY+x][SHIFTXY+y]==1 ) {
							 if( (s=Tuning_Image[SHIFTXY+x+1][SHIFTXY+y  ])<0 ) { Segmented_Image[SHIFTXY+x][SHIFTXY+y]=s; n++; }
						else if( (s=Tuning_Image[SHIFTXY+x-1][SHIFTXY+y  ])<0 ) { Segmented_Image[SHIFTXY+x][SHIFTXY+y]=s; n++; }
						else if( (s=Tuning_Image[SHIFTXY+x  ][SHIFTXY+y+1])<0 ) { Segmented_Image[SHIFTXY+x][SHIFTXY+y]=s; n++; }
						else if( (s=Tuning_Image[SHIFTXY+x  ][SHIFTXY+y-1])<0 ) { Segmented_Image[SHIFTXY+x][SHIFTXY+y]=s; n++; }
						else Segmented_Image[SHIFTXY+x][SHIFTXY+y] = Tuning_Image[SHIFTXY+x][SHIFTXY+y];
					} else   Segmented_Image[SHIFTXY+x][SHIFTXY+y] = Tuning_Image[SHIFTXY+x][SHIFTXY+y];
		} while( n>0 ); 
	}

	// Numbering REMAINING segmented areas with increasingly NEGATIVE numbers  
	for(x=0; x<Bimage->xsize; x++) 		
		for(y=0; y<Bimage->ysize; y++)
			if( (s=Segmented_Image[SHIFTXY+x][SHIFTXY+y]) == 1 ){						// Pixel identified as dark region not overlapping with any previous cell
				a = fill_segmented(s, -counter, x, y);									// Filling cell with NEGATIVE number
				counter++;
			}

/*	if( Write_4_CU_img ) {
		copy_bitmap(Rimage, Oimage);								// Copying original image into output image 
		convert2color(Oimage);										// To allow showing different cells with different colors 
		for( y=0; y<Bimage->ysize; y++ )
			for( x=0; x<Bimage->xsize; x++ ) 
				if( (s=Segmented_Image[SHIFTXY+x][SHIFTXY+y])<0 && (Segmented_Image[SHIFTXY+x+1][SHIFTXY+y]!=s || Segmented_Image[SHIFTXY+x-1][SHIFTXY+y]!=s || Segmented_Image[SHIFTXY+x][SHIFTXY+y+1]!=s || Segmented_Image[SHIFTXY+x][SHIFTXY+y-1]!=s ) ) {
					cellnumber2color(s,color);			// The colorcode sets channel by bit number 1=blue 2=green 4=red. Can be combined to give shades of cyan=3, magenta=5, yellow=6. 
					for( i=0; i<3; i++) Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][i] = color[i]; 			
				}
		if( Add_Scalebar ) draw_scalebar(); 
		sprintf(filename, "%s/%s_04_3_%04d_New-regions-filled.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*level);
		out_bitmap(Oimage, filename);	
	} //*/

/*	// For images after the zeroth, making the copy of the Previous_Segmented_Image[][] here in order to keep small and large candidate particles which are scrapped below 
	if( level > 0 ) 
		for(x=0; x<Bimage->xsize; x++) 		
			for(y=0; y<Bimage->ysize; y++)
				Previous_Segmented_Image[SHIFTXY+x][SHIFTXY+y] = Segmented_Image[SHIFTXY+x][SHIFTXY+y];		//*/
	
	if( Perform_Tuning )								// Making a copy of the segmented cell image before scrapping by size below to allow displaying too small or large regions in tuning image
		for(x=0; x<Bimage->xsize; x++) 		
			for(y=0; y<Bimage->ysize; y++)
				Tuning_Image[SHIFTXY+x][SHIFTXY+y] = Segmented_Image[SHIFTXY+x][SHIFTXY+y];			// Tuning_Image[][] thus now has all segmented regions coded as negative numbers

	// Identifying individual cells and numbering them in the Cells[][] array and Segmented_Image[][] with POSITIVE numbers
	Nof__Cells[level]=0;																							// NB: Nof__Cells[] counts from 1 (not zero) in order to match the fill number in Segmented_Image[][] where zero is busy as background
	for(x=0; x<Bimage->xsize; x++) 		
		for(y=0; y<Bimage->ysize; y++)
			if( (s=Segmented_Image[SHIFTXY+x][SHIFTXY+y]) < 0 ){													// Pixel belongs to region marked with NEGATIVE number above  
				i = ++Nof__Cells[level];
				if( i == MAXCELLS-1 ) { sprintf(message, "The maximum number of cells per image (%d) was reached at image number %d.\nReduce the number of images or increase MAXCELLS and recompile celltraxx.c.\n", MAXCELLS, First_Image_Number+Image_Nr_Increment*level);  write_error_message(ErrorMessageFile, ERROR, message); }
				Cells[level][Nof__Cells[level]].cm_x_first = x;														// Keeping this point as a safe place to start filling when matching cells below P
				Cells[level][Nof__Cells[level]].cm_y_first = y;
				Cells[level][Nof__Cells[level]].cellnumber = (int)(Nof__Cells[level]);								// Temporarily assigning a cell number to this cell to allow outputting image ...05_Identified.bmp with different colors 
				a = Cells[level][Nof__Cells[level]].area_px = fill_segmented_box(s, Nof__Cells[level], x, y, &Cells[level][Nof__Cells[level]].min_x, &Cells[level][Nof__Cells[level]].max_x, &Cells[level][Nof__Cells[level]].min_y, &Cells[level][Nof__Cells[level]].max_y);			// Filling cell with POSITIVE identification number and storing number of pixels in .area_px
				hist_index = (int)( Nof_Bins*2*sqrt(a/3.141592)/Max_Diameter_px );									// Equivalent circular diameter index when  Diameter_Histogram[][] goes from zero to Max_Diameters
				box_width  = Cells[level][Nof__Cells[level]].max_x - Cells[level][Nof__Cells[level]].min_x;		 
				box_height = Cells[level][Nof__Cells[level]].max_y - Cells[level][Nof__Cells[level]].min_y;		 	 
				box_midpoint_y = (int)(0.5*( Cells[level][Nof__Cells[level]].max_y + Cells[level][Nof__Cells[level]].min_y ));		
//				if( a<Min_Cell_Pixels || box_width<Min_Diameter_px || box_height<Min_Diameter_px ) { 
				if( a<Min_Cell_Pixels ) { 
					fill_segmented(Nof__Cells[level], 0, x, y);														// Filling regions with background color and correcting counter if region is too small
					Nof__Cells[level]--;	
					if( Perform_Tuning ) fill_tuning(s, TOO_SMALL_CELL, x, y);										// Coding too small regions as +TOO_SMALL_CELL in Tuning_Image[][]
//				} else if( a>Max_Cell_Pixels || box_width>Max_Diameter_px || box_height>Max_Diameter_px ) { 
				} else if( a>Max_Cell_Pixels ) { 
					fill_segmented(Nof__Cells[level], 0, x, y);														// Filling regions with background color and correcting counter if region is too large 
					Nof__Cells[level]--;								
					if( Perform_Tuning ) fill_tuning(s, TOO_LARGE_CELL, x, y);										// Coding too large regions as +TOO_LARGE_CELL in Tuning_Image[][]
				} else if( Wound_Healing_Mode && level>0 && ( box_midpoint_y>WH_Box_Upper_Limit_px || box_midpoint_y<WH_Box_Lower_Limit_px ) ) {																	// In wound healing mode, removing all cells that are too far from the wound front(s)
					fill_segmented(Nof__Cells[level], 0, x, y);														// Filling regions with background color and correcting counter if region is too large 
					Nof__Cells[level]--;								
					if( Perform_Tuning ) fill_tuning(s, TOO_FAR_FROM_WH_EDGE, x, y);								// Coding regions too far from the wound healing front as +TOO_FAR_FROM_WH_EDGE in Tuning_Image[][]		
				} else Diameter_Histogram[1][hist_index]++;															// Counting all accepted  cells in second column of Diameter_Histogram[1][]
				Diameter_Histogram[0][hist_index]++; 																// Counting all candidate cells in first  column of Diameter_Histogram[0][]
				if( box_height<=0 || box_width<=0 ) {
					Cells[level][Nof__Cells[level]].aspect_ratio  = 0;		// Just coding as zero to avoid 1.#INF0 in PosFile
					Cells[level][Nof__Cells[level]].area_fraction = 0; 
				} else {
					Cells[level][Nof__Cells[level]].aspect_ratio  = 1.0*box_height/box_width; 
					Cells[level][Nof__Cells[level]].area_fraction = 4.0*a/(Pi*box_height*box_width); 
				}
			} 

	// Copying to Previous_Segmented_Image[][] here after scrapping small and large candidate particles above 
	for(x=0; x<Bimage->xsize; x++) 		
		for(y=0; y<Bimage->ysize; y++)
			Previous_Segmented_Image[SHIFTXY+x][SHIFTXY+y] = Segmented_Image[SHIFTXY+x][SHIFTXY+y];		
	
	// Writing IDENTIFIED image to file for checking. Drawing only cell contours for more easy comparison 
	if( Write_5_ID_img || Make__5_ID_vid ) {
		copy_bitmap(Rimage, Oimage);								// Copying original image into output image 
		convert2color(Oimage);										// To allow showing different cells with different colors 
		for( y=0; y<Bimage->ysize; y++ )
			for( x=0; x<Bimage->xsize; x++ ) 
				if( Segmented_Image[SHIFTXY+x][SHIFTXY+y] == CUT ) {
					Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0] = Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1] = (unsigned char)0; Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2] = (unsigned char)255;  			
				} else if( (s=Segmented_Image[SHIFTXY+x][SHIFTXY+y])>0 && (Segmented_Image[SHIFTXY+x+1][SHIFTXY+y]!=s || Segmented_Image[SHIFTXY+x-1][SHIFTXY+y]!=s || Segmented_Image[SHIFTXY+x][SHIFTXY+y+1]!=s || Segmented_Image[SHIFTXY+x][SHIFTXY+y-1]!=s ) ) {
					cellnumber2monocolor(s,color,2);	// The colorcode sets channel by bit number 1=blue 2=green 4=red. Can be combined to give shades of cyan=3, magenta=5, yellow=6. 
					for( i=0; i<3; i++) Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][i] = color[i]; 			
				}
		if( Add_Scalebar ) draw_scalebar(); 
		sprintf(filename, "%s/%s_05_%04d_Identified.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*level);
		out_bitmap(Oimage, filename);
	}

	// Option to write IDENTIFIED MASK image to file for mask segmentation and tracking in other software. White cells (255) on black background (0).
/*	if( TRUE ) {
		copy_bitmap(Rimage, Oimage);								// Copying original image into output image 
		for( y=0; y<Bimage->ysize; y++ )
			for( x=0; x<Bimage->xsize; x++ )				// Including only pixels that are surrounded by the same cell region pixels on all sides to split touching cells 
				if( (s=Segmented_Image[SHIFTXY+x][SHIFTXY+y])>0 && Segmented_Image[SHIFTXY+x+1][SHIFTXY+y]==s && Segmented_Image[SHIFTXY+x-1][SHIFTXY+y]==s && Segmented_Image[SHIFTXY+x][SHIFTXY+y+1]==s && Segmented_Image[SHIFTXY+x][SHIFTXY+y-1]==s  )	
						Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0] = Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1] = Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2] = (unsigned char)255;  			// Have to set all three channels even if not color image since channel averaging is used when saving 8-bit images 		
				else	Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0] = Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1] = Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2] = (unsigned char)0; 			
		sprintf(filename, "%s/%s__5_%04d_Identified_Binary.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*level);
		out_bitmap(Oimage, filename);
	} //*/
	
	// Calculating the center of mass coordinates, in parallel for all cells 
	for( y=0; y<Bimage->ysize; y++ )
		for( x=0; x<Bimage->xsize; x++ ) {
			Cells[level][ Segmented_Image[SHIFTXY+x][SHIFTXY+y] ].dx += x;
			Cells[level][ Segmented_Image[SHIFTXY+x][SHIFTXY+y] ].dy += y;
		}
	for( c=1; c<=Nof__Cells[level]; c++ ) {
		Cells[level][c].dx /= Cells[level][c].area_px;			// Dividing the summed moments by the area to get the center of mass coordinate 
		Cells[level][c].dy /= Cells[level][c].area_px;
		Cells[level][c].cx  = (int)(Cells[level][c].dx+0.5);	// Rounding off the integer version of the center of mass coordinates
		Cells[level][c].cy  = (int)(Cells[level][c].dy+0.5);		
	}
}
	

// Splitting convex regions by finding all edge points and their mutual distances, then drawing cutting lines to erase segmented pixels between narrowly spaced spaced points where at least one side is convex 

char cut_concave_region(int start_x, int start_y, int level)
{
	char found;	
	short ss; 
	int e, f, i, j, k, s, w, x, y, xmin, xmax, ymin, ymax, xend, yend, xmid, ymid, nof_edge_points, max_nof_edge_points;
	double ri, rj; 
	double step_x, step_y, xd, yd, dist; 
//	char filename[200];	FILE *csvfile;	// Option to write selected Edge_Point_Matrix[][]s to file for checking shape, edge point numbering, etc.

	ss = (short)Nof__Cells[level];			// Positive number filled into current region to find region size
	xmin = Cells[level][ss].min_x-1;		// Shorthand + Expanding borders by 1 pixel on each side so no edge points are on border of Edge_Point_Matrix[][]
	xmax = Cells[level][ss].max_x+1;		 
	ymin = Cells[level][ss].min_y-1;
	ymax = Cells[level][ss].max_y+1;

	start_x = start_x-xmin;					// Changing coordinates to match Edge_Point_Matrix[][]
	start_y = start_y-ymin;
	xend = xmax-xmin;
	yend = ymax-ymin;

//	printf("Running cut_concave_region() with x=%5d y=%5d at level=%4d with Nof__Cells = %4d   xend = %4d   yend = %4d \n", start_x, start_y, level, ss, xend, yend); 

	if( xend > MAXX ) return(FALSE);	// Should never happen, but just checking anyway
	if( yend > MAXX ) return(FALSE); 

	// Filling the local matrix Edge_Point_Matrix[][] with -1 at edge points and 0 at other points for this region
	for( x=xmin; x<=xmax; x++ ) 
		for( y=ymin; y<=ymax; y++ ) 
			if( Segmented_Image[SHIFTXY+x][SHIFTXY+y]!=ss && (Segmented_Image[SHIFTXY+x+1][SHIFTXY+y]==ss || Segmented_Image[SHIFTXY+x-1][SHIFTXY+y]==ss || Segmented_Image[SHIFTXY+x][SHIFTXY+y+1]==ss || Segmented_Image[SHIFTXY+x][SHIFTXY+y-1]==ss ) ) 
				 Edge_Point_Matrix[x-xmin][y-ymin] = -1;	// Coding edge points as negative one
			else Edge_Point_Matrix[x-xmin][y-ymin] =  0;

	// Searching all around the region and storing edge points in coordinate arrays Edge_Point_X_px[] and Edge_Point_Y_px[] 
	max_nof_edge_points = CUT_MAX_EDGE_POINTS - 2*CUT_MARGIN - 2*CIRCLE_FIT_POINT_DIST; 
	nof_edge_points = 0;
	Edge_Point_X_px[0] = x = start_x-1;		// Special treatment for 0th point, stepping one pixel left to the column that does not contain any segmented pixels 
	Edge_Point_Y_px[0] = y = start_y;

	do {
		found = FALSE; 
		e = 1;	// Edge length in spiral
		s = 1;	// Sign for movement in spiral 
		// Searching outward in "square shaped spiral" for the next edge point, beginning upwards and to the right
		do{ 
			f = 0; while( !found && f<e ) { f++; y += s; if( x>=0 && x<=xend && y>=0 && y<=yend && Edge_Point_Matrix[x][y]==-1 ) found=TRUE; }
			f = 0; while( !found && f<e ) { f++; x += s; if( x>=0 && x<=xend && y>=0 && y<=yend && Edge_Point_Matrix[x][y]==-1 ) found=TRUE; }
			s *= -1; 
		} while( !found && ++e<MAX_SEARCH_SPIRAL_EDGE );	// Searching 5 pixels to each side (when MAX_SEARCH_SPIRAL_EDGE = 9). Should be enough in most cases. 

		if( !found ) return(FALSE); // { printf("Error in cut_conca ve_region(): Could not find the next edge point even if the search spiral \nreached an edge length of %d pixelsfor the region beginning at (%4d,%4d).\nPlease increases MAX_SEARCH_SPIRAL_EDGE and recompile.\n", e, start_x+xmin, start_y+ymin); exit(TRUE); }
		else if( x>=0 && x<=xend && y>=0 && y<=yend && x<MAXX && y<MAXY ) {		// Double checking that point is valid and that array limits will not be overwritten.
			if( ++nof_edge_points == max_nof_edge_points ) { printf("\nAborting cut_concave_region() since maxumum allowed number of edge points (%d) was reached for region starting at pixel [%d,%d]\nConsider increasing CUT_MAX_EDGE_POINTS and recompiling celltraxx.c if the segmented image looks good.\n", max_nof_edge_points, start_x+xmin, start_y+ymin); return(FALSE); }  
			Edge_Point_X_px[nof_edge_points] = x;	
			Edge_Point_Y_px[nof_edge_points] = y;
			Edge_Point_Matrix[x][y] = nof_edge_points;							// Giving point a new number so not counted again
		}   
	} while( x!=Edge_Point_X_px[0] || y!=Edge_Point_Y_px[0] );

	// Extending point list to get overlap in below matrix 
	for( i=1; i<2*CIRCLE_FIT_POINT_DIST; i++ ) {
		Edge_Point_X_px[nof_edge_points+i] = Edge_Point_X_px[i];
		Edge_Point_Y_px[nof_edge_points+i] = Edge_Point_Y_px[i];
	}

/*	if( nof_edge_points > 300 ) {			// Option to write selected Edge_Point_Matrix[][]s to file for checking shape, edge point numbering, etc.
		printf("\n nof_edge_points = %4d at [%4d,%4d]   xmin=%d  xmax=%d  ymin=%d  ymax=%d  ", nof_edge_points, start_x+xmin, start_y+ymin, xmin, xmax, ymin, ymax); 
		sprintf(filename, "%s/%s_05_%04d_csvfile_x=%04d_y=%04d.csv", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*level, start_x+xmin, start_y+ymin); 
		csvfile = fopen(filename, "w"); 
		for( y=ymax; y>=ymin; y-- ) {
			for( x=xmin; x<=xmax; x++ ) fprintf(csvfile,"%d,", Edge_Point_Matrix[x-xmin][y-ymin] );
			fprintf(csvfile,"\n"); 
		}
		fprintf(csvfile,"\n");
		fclose(csvfile);
	} //*/

	// Generating upper triangle in matrix with distances between all edge points 
	for( i=0; i<nof_edge_points+2*CIRCLE_FIT_POINT_DIST; i++ ) 
		for( j=i; j<nof_edge_points+2*CIRCLE_FIT_POINT_DIST; j++ )   
			Edge_Point_Distances[i][j] = sqrt( 1.0*(Edge_Point_X_px[j]-Edge_Point_X_px[i])*(Edge_Point_X_px[j]-Edge_Point_X_px[i]) + (Edge_Point_Y_px[j]-Edge_Point_Y_px[i])*(Edge_Point_Y_px[j]-Edge_Point_Y_px[i]) );	

	// Finding local minima and drawing lines in Segmented_Image[][] to separate regions between the minima.	NB: To avoid errors below, never use indeces in Edge_Point_Distances[][] that are more than ± CUT_MARGIN away from i or j.
	for( i=CUT_MARGIN; i<nof_edge_points+CUT_MARGIN; i++ ) {
		for( j=i+CUT_MARGIN; j<nof_edge_points+CUT_MARGIN; j++ ) {								// Starting search CUT_MARGIN points past i since no reason to begin any nearer. NB: To avoid errors below, never use indeces in Edge_Point_Distances[][] that are more than ± CUT_MARGIN away from i or j.
			dist=Edge_Point_Distances[i][j];
			found = TRUE;
			for( k=-CUT_MARGIN; k<=CUT_MARGIN; k++ ) if( Edge_Point_Distances[i+k][j  ] < dist ) found = FALSE;
			for( k=-CUT_MARGIN; k<=CUT_MARGIN; k++ ) if( Edge_Point_Distances[i  ][j+k] < dist ) found = FALSE;
			if( i > nof_edge_points ) e = i-nof_edge_points; else e = i;
			if( j > nof_edge_points ) f = j-nof_edge_points; else f = j;
			if( found && mini(abs(f-e), abs(e-f+nof_edge_points)) > 0.7*Pi*Cutting_Diameter_px ) {		// Last test finds smallest number of edge points between i and j, and avoids cutting if the cut protrusion would have a smaller arch length than 80 % of the typical cutting cell circumference. 

				ri = radius_of_curvature_for_fitted_circle(i); 
				rj = radius_of_curvature_for_fitted_circle(j); 

				xmid = (int)( 0.5*(Edge_Point_X_px[i]+Edge_Point_X_px[j]) + 0.5 );						// Midpoint between the two narrow point, used below to make sure this point is inside the cell region so the cut will actually divide the region
				ymid = (int)( 0.5*(Edge_Point_Y_px[i]+Edge_Point_Y_px[j]) + 0.5 );

//				if( (dist<5 || ri+rj<20) && Segmented_Image[SHIFTXY+xmin+xmid][SHIFTXY+ymin+ymid]==ss ) {								// Demanding a short separation or a small radius of curvature on at least one of the sides  
//				if( dist*(ri+rj)<500 && Segmented_Image[SHIFTXY+xmin+xmid][SHIFTXY+ymin+ymid]==ss ) {								// Demanding a short separation or a small radius of curvature on at least one of the sides  
//				if( dist*(ri+rj)<500 || ri<8 || rj<8 ) {								// Demanding a short separation or a small radius of curvature on at least one of the sides  
//				if( ri<BIG-1 && rj<BIG-1 && dist*mind(ri,rj)<100 && (ri<30 || rj<30) && Segmented_Image[SHIFTXY+xmin+xmid][SHIFTXY+ymin+ymid]==ss ) {								// Demanding a short separation or a small radius of curvature on at least one of the sides  
//				if(  dist*dist*mind(ri,rj)<2000 && ri<BIG-1 && rj<BIG-1 && Segmented_Image[SHIFTXY+xmin+xmid][SHIFTXY+ymin+ymid]==ss ) {								// Demanding a short separation or a small radius of curvature on at least one of the sides  
//				if( ri<BIG-1 && rj<BIG-1 && dist*mind(ri,rj)<Cutting_Diameter_px && Segmented_Image[SHIFTXY+xmin+xmid][SHIFTXY+ymin+ymid]==ss ) {			// Demanding both points are concave AND a short separation or a small radius of curvature on at least one of the sides AND that the midpoint of the cut is within the region
//				if( (ri<BIG-1 || rj<BIG-1) && dist*mind(ri,rj)<Cutting_Diameter_px && Segmented_Image[SHIFTXY+xmin+xmid][SHIFTXY+ymin+ymid]==ss ) {			// Demanding both points are concave AND a short separation or a small radius of curvature on at least one of the sides AND that the midpoint of the cut is within the region
//				if( ( ri>0 && ri<99 && ( rj>0 || rj<-99) ) && dist*mind(abs(ri),abs(rj))<Cutting_Diameter_px && Segmented_Image[SHIFTXY+xmin+xmid][SHIFTXY+ymin+ymid]==ss ) {			// Demanding both points are concave AND a short separation or a small radius of curvature on at least one of the sides AND that the midpoint of the cut is within the region
//				if( ( ri>0 && ri<99 && ( rj>0 || rj<-99) ) && dist*mind(abs(ri),abs(rj))<Cutting_Diameter_px && Segmented_Image[SHIFTXY+xmin+xmid][SHIFTXY+ymin+ymid]==ss ) {			// Demanding both points are concave AND a short separation or a small radius of curvature on at least one of the sides AND that the midpoint of the cut is within the region
//				if( ( ri>0 && ri<999 && ( rj>0 || rj<-99) ) && dist*mind(abs(ri),abs(rj))<Cutting_Diameter_px && Segmented_Image[SHIFTXY+xmin+xmid][SHIFTXY+ymin+ymid]==ss ) {			// Demanding both points are concave AND a short separation or a small radius of curvature on at least one of the sides AND that the midpoint of the cut is within the region
//				if( ri>0 && rj>0 && dist*mind(abs(ri),abs(rj))<Cutting_Diameter_px && Segmented_Image[SHIFTXY+xmin+xmid][SHIFTXY+ymin+ymid]==ss ) {			// Demanding both points are concave AND a short separation or a small radius of curvature on at least one of the sides AND that the midpoint of the cut is within the region
				if( (( ri>0 && rj>0 && ri+rj<2.5*Cutting_Diameter_px && dist<Cutting_Diameter_px ) || ( ri>0 && rj>0 && dist<0.333*Cutting_Diameter_px )) && Segmented_Image[SHIFTXY+xmin+xmid][SHIFTXY+ymin+ymid]==ss ) {			// Demanding both points are concave AND a short separation or a small radius of curvature on at least one of the sides AND that the midpoint of the cut is within the region
//				if( ri>0 && rj>0 ) {			// Demanding both points are concave AND a short separation or a small radius of curvature on at least one of the sides AND that the midpoint of the cut is within the region
//				if( ri>0 && rj!=0 && dist>0 && 60000*(1/ri + 1/rj)/(dist*dist) > Cutting_Diameter_px && Segmented_Image[SHIFTXY+xmin+xmid][SHIFTXY+ymin+ymid]==ss ) {			// Demanding both points are concave AND a short separation or a small radius of curvature on at least one of the sides AND that the midpoint of the cut is within the region
					step_x = 1.0*(Edge_Point_X_px[j]-Edge_Point_X_px[i])/dist; 
					step_y = 1.0*(Edge_Point_Y_px[j]-Edge_Point_Y_px[i])/dist; 
					for( w = -1; w <= 1; w++ ) {									// Erasing lines that are 2 pixels wide 
						xd = Edge_Point_X_px[i] + 0.6*w*step_y -1*step_x;			// The exact, fractional coordinates of the pixels. Starting 1 pixels beoyond the narrowest point to get a clean cut.
						yd = Edge_Point_Y_px[i] - 0.6*w*step_x -1*step_y;				
						for( k=0; k<=dist+2; k++ ) {								// Ending 3 pixels beoyond the narrowest point to get a clean cut.	//*/
							x = (int)(xd + 0.5);									// The rounded off integer pixel coordinates to use as array indeces
							y = (int)(yd + 0.5);
							Segmented_Image[SHIFTXY+xmin+x][SHIFTXY+ymin+y] = Matched_Image[SHIFTXY+xmin+x][SHIFTXY+ymin+y] = CUT;		// Using Matched_Image[][] as temporary storage for cut lines such that they can be displayed in the tuning image. Needed since the cut lines and other FIT markings are deleted from Segmented_Image[][] for further processing.
							xd += step_x; 
							yd += step_y; 
						}
					}
				} 
				if( ri>0 ) Segmented_Image[SHIFTXY+xmin+Edge_Point_X_px[i]][SHIFTXY+ymin+Edge_Point_Y_px[i]] = FIT; else if( ri<0 ) Segmented_Image[SHIFTXY+xmin+Edge_Point_X_px[i]][SHIFTXY+ymin+Edge_Point_Y_px[i]] = FIT2; else Segmented_Image[SHIFTXY+xmin+Edge_Point_X_px[i]][SHIFTXY+ymin+Edge_Point_Y_px[i]] = FIT3; 
				if( rj>0 ) Segmented_Image[SHIFTXY+xmin+Edge_Point_X_px[j]][SHIFTXY+ymin+Edge_Point_Y_px[j]] = FIT; else if( rj<0 ) Segmented_Image[SHIFTXY+xmin+Edge_Point_X_px[j]][SHIFTXY+ymin+Edge_Point_Y_px[j]] = FIT2; else Segmented_Image[SHIFTXY+xmin+Edge_Point_X_px[j]][SHIFTXY+ymin+Edge_Point_Y_px[j]] = FIT3;
//				printf("Finished cutting between edge points %d and %d \n", i, j); 
			}
		}
	}	//*/
	return(TRUE); 
}


// Calculating the radius of curvature of a circle fitted to point number i on the edge of a region and its neighbours ±CIRCLE_FIT_POINT_DISTs away

double radius_of_curvature_for_fitted_circle(int i) 
{
	double x1, y1, x2, y2, x3, y3, x4, y4, d1, d2, a, b, r, z; 

	x1 = Edge_Point_X_px[i-CIRCLE_FIT_POINT_DIST];						 
	y1 = Edge_Point_Y_px[i-CIRCLE_FIT_POINT_DIST]; 
	x2 = Edge_Point_X_px[i]; 
	y2 = Edge_Point_Y_px[i]; 
	x3 = Edge_Point_X_px[i+CIRCLE_FIT_POINT_DIST]; 
	y3 = Edge_Point_Y_px[i+CIRCLE_FIT_POINT_DIST]; 

	z = (x3-x2)*(y1-y2)-(y3-y2)*(x1-x2);													// The z component of the cross product 23 x 21 which will be positive in concave regions, negative in convex regions and zero if the three vectors are co-linear
	if( z == 0 ) return(-999);																// If three points in a straight line, just returning a large, convex radius of curvature

	d2 = 2*(x2-x1);		
	if( d2==0 ) { x4=x1; y4=y1; x1=x3; y1=y3; x3=x4; y3=y4; d2 = 2*(x2-x1); }				// If point 1 and 2 happen to have the same x coordinate, swapping points 1 and 3 and recalculating d2 since point 2 and 3 should not have the same x coordinate since z is not zero here and the three points are not co-linear 

	d1 = 2*((x3-x1)*(y1-y2)-(x1-x2)*(y3-y1)); 												// Twice the z component of the cross product 13 x 21 which will also be zero if the three vectors are co-linear 
	if( d1==0 || d2==0 ) return(0); 														// Should not happen after the above tests and swapping, but just double checking to guard against division by zero error  
	else {
		b = ((x2-x1)*(x3*x3+y3*y3-x1*x1-y1*y1)-(x3-x1)*(x2*x2+y2*y2-x1*x1-y1*y1))/d1;		// Center y coordinate of circle
		a = ((x2*x2+y2*y2-x1*x1-y1*y1)+2*b*(y1-y2))/d2;										// Center x coordinate of circle
		r = sqrt((x1-a)*(x1-a)+(y1-b)*(y1-b));
		if( z < 0  ) return(-r);															// Returning negative radius if convex point 
		else         return( r); 
	}
}


// Auxiliary function used in cut_concave_region();

double mind(double a, double b)
{
	if( a < b ) return(a);
	else        return(b);
}


// For the first image in the first video, write to disk the result of the current analysis with blue outlines for accepted, yellow for too small cells and red for too large cells. Image is used by ImageJ macro for interactive tuning.

void write_interactively_tuned_image(void) 
{
	char filename[STRING_LENGTH]; 
	unsigned char color[3]; 
	short s; 
	int i, x, y, y_shift; 

	copy_bitmap(Rimage, Oimage);													// Copying original image into output image 
	convert2color(Oimage);															// To allow showing different cells with different color outlines 

	for( y=0; y<Bimage->ysize; y++ )												// Tuning_Image[][] has too small/large regions coded as +TOO_SMALL/LARGE_CELL and valid regions coded as negative numbers 
		for( x=0; x<Bimage->xsize; x++ ) {
			s = Tuning_Image[SHIFTXY+x][SHIFTXY+y];
 			if(      s < 0                     && (Tuning_Image[SHIFTXY+x+1][SHIFTXY+y]!=s || Tuning_Image[SHIFTXY+x-1][SHIFTXY+y]!=s || Tuning_Image[SHIFTXY+x][SHIFTXY+y+1]!=s || Tuning_Image[SHIFTXY+x][SHIFTXY+y-1]!=s ) ) { cellnumber2monocolor(-s,color,2); for( i=0; i<3; i++)  Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][i] = color[i]; }  // Trick to get different shades of green for accepted cells to better see if neighbours are separate or not 	
			else if( s == TOO_SMALL_CELL       && (Tuning_Image[SHIFTXY+x+1][SHIFTXY+y]!=s || Tuning_Image[SHIFTXY+x-1][SHIFTXY+y]!=s || Tuning_Image[SHIFTXY+x][SHIFTXY+y+1]!=s || Tuning_Image[SHIFTXY+x][SHIFTXY+y-1]!=s ) ) { Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0] = (unsigned char)(  0); Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1] = (unsigned char)(255); Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2] = (unsigned char)(255); }  		
			else if( s == TOO_LARGE_CELL       && (Tuning_Image[SHIFTXY+x+1][SHIFTXY+y]!=s || Tuning_Image[SHIFTXY+x-1][SHIFTXY+y]!=s || Tuning_Image[SHIFTXY+x][SHIFTXY+y+1]!=s || Tuning_Image[SHIFTXY+x][SHIFTXY+y-1]!=s ) ) { Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0] = (unsigned char)(  0); Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1] = (unsigned char)(  0); Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2] = (unsigned char)(255); }  		
			else if( s == TOO_FAR_FROM_WH_EDGE && (Tuning_Image[SHIFTXY+x+1][SHIFTXY+y]!=s || Tuning_Image[SHIFTXY+x-1][SHIFTXY+y]!=s || Tuning_Image[SHIFTXY+x][SHIFTXY+y+1]!=s || Tuning_Image[SHIFTXY+x][SHIFTXY+y-1]!=s ) ) { Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0] = (unsigned char)(  0); Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1] = (unsigned char)(  0); Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2] = (unsigned char)(  0); }  		
			if( Matched_Image[SHIFTXY+x][SHIFTXY+y] == CUT )	/* Cut lines are temporarily coded as CUT in Matched_Image[][] to allow display here */																			{ Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0] = (unsigned char)(255); Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][1] = (unsigned char)(  0); Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][2] = (unsigned char)(  0); } 		
		}

	sprintf(filename, "C:/celltraxx_system/celltraxx_tuning_image_inset.bmp");						// Image made in Paint with font Arial, size 24 for black text and 21 for colored text
	read_bitmap(Dimage, filename);	
	Oimage->ysize += Dimage->ysize;																	// Expanding the output image with the size of the bottom part image to allow pasting it in below
	y_shift = Oimage->ysize - Dimage->ysize;														

	for( y=0; y<Dimage->ysize; y++ )						
		for( x=0; x<Dimage->xsize; x++ )															// ...then pasting inset image at the top left of expanded image for output to ImageJ macro
			for( i=0; i<3; i++ ) Oimage->image[x][y+y_shift][i] = Dimage->image[x][y][i];

	for( i=0; i<3; i++ ) color[i] = Dimage->image[Dimage->xsize-1][0][i];							// Copying color of bottom right pixel in inset image 
	for( y=0; y<Dimage->ysize; y++ )						
		for( x=Dimage->xsize; x<Oimage->xsize; x++ )												// Just filling upper right corner of expanded image with the same color (to avoid black rectangle) 
			for( i=0; i<3; i++ ) Oimage->image[x][y+y_shift][i] = color[i];

	if( Add_Scalebar ) draw_scalebar(); 
	sprintf(filename, "C:/celltraxx_system/tuning_temp/celltraxx_tuning_image.bmp");
	out_bitmap(Oimage, filename);	
/*
	// Writing updated info back to the tuning image file as a way to remember the segmentation limit if the same image is analysed again and if tuning the middle and last image of the same video.
	sprintf(filename, "C:/celltraxx_system/celltraxx_tuning_info.txt"); 
	if( (TuneFile=fopen(filename, "w"))==NULL ){ sprintf(message, "Can't open file celltraxx_tuning_info.txt for writing.\nPlease check if the file is open in another application.\nExiting.\n");  write_error_message(ErrorMessageFile, ERROR, message); }
	fprintf(TuneFile, "First image grey level limit %5.3f\n", Segmentation_Grey_Limit_From_First_Image);
	fprintf(TuneFile, "First image path and name    %s\n",    PreviousFirstImageName);
	fclose(TuneFile);	//*/
}


// Finding the most likely match of cells in the current level and the previous level by matching first the cells that have moved the least 

void match_cells_with_previous_level(int level) 
{
	char cells_matched, filename[STRING_LENGTH], message[STRING_LENGTH]; 
	unsigned char color[3];
	short s; 
	int c, d, i, j, i_1, j_1, j_2, x, y, int_radius; 
	int x_1, y_1, x_2, y_2, m, m_1, m_2;	
	int delta_x, delta_y, delta_max, lev; 
	float shift, smallest_shift_1, dot_radius; 
	float smallest_shift_2=0; 
	double core=0, alpha=0; 

	// Special treatment for zeroth level which does not have a previous level to estimate cell movements from 
	if( level == 0 ){
		for( c=1; c<=Nof__Cells[0]; c++ ) {					
			Cells[0][c].cellnumber = c;												// Assigning each physical cell identified in the first image a unique number, which will then be assigned also to Cells[][] at other levels (which will typically have different last index value) 
			Matched_Cells[0][c] = c;												// Also populating the first level of the Matched_Cells[][] overview matrix with the Nof__Cells[1] cells found in the first image
			fill_segmented(c, -c, Cells[0][c].cm_x_first, Cells[0][c].cm_y_first);	// Marking all cells in the first level as "matched" by filling it with the NEGATIVE of its cell number, to be compatible with same marking for cells in subsequent levels 

/*			// Temporarily drawing bounding boxes around all matched cells  
			for( x=Cells[level][c].min_x; x<=Cells[level][c].max_x; x++ ) Segmented_Image[SHIFTXY+x][SHIFTXY+Cells[level][c].min_y] = Segmented_Image[SHIFTXY+x][SHIFTXY+Cells[level][c].max_y] = -c;						
			for( y=Cells[level][c].min_y; y<=Cells[level][c].max_y; y++ ) Segmented_Image[SHIFTXY+Cells[level][c].min_x][SHIFTXY+y] = Segmented_Image[SHIFTXY+Cells[level][c].max_x][SHIFTXY+y] = -c;			//*/ 			
		}
		Nof_Matched_Cells = Nof_Identified_Cells[0] = Nof__Cells[0]; 
		printf(         "%5d     %5d \n", Nof__Cells[0], 0);
		fprintf(CMDWindowFile, "%d,%d\n", Nof__Cells[0], 0);

	} else {	// From level 1 and higher levels 

//		printf("\nNof__Cells[level-1] = %4d    Nof__Cells[level] = %4d \n", Nof__Cells[level-1], Nof__Cells[level]); 
		// Populating the Matching_Matrix[][] with the distances between all cells in the previous level and this level, in order to match the cells with the smallest distance first
		for( i=1; i<=Nof__Cells[level-1]; i++ )
			for( j=1; j<=Nof__Cells[level]; j++ )
				Matching_Matrix[i][j] = (float)sqrt(     (Cells[level][j].dx-Cells[level-1][i].dx)*(Cells[level][j].dx-Cells[level-1][i].dx) + (Cells[level][j].dy-Cells[level-1][i].dy)*(Cells[level][j].dy-Cells[level-1][i].dy) );
				// USING ONLY DISTANCE AS MATCHING CRITERION FOR NOW. CONSIDER ADDING AREA, ROUNDNESS, TILT ANGLE ETC BY IMPLEMENTING ELLIPSE FITTING LATER 

		for( i=0; i<=Nof__Cells[level-1]; i++ ) Matching_Matrix[i][0]        = (float)i;		// Filling zeroth  row   with cell numbers  
		for( j=1; j<=Nof__Cells[level  ]; j++ ) Matching_Matrix[0][j]        = (float)j;		// Filling zeroth column with cell numbers  
		for( j=1; j<=Nof__Cells[level  ]; j++ ) Matching_Matrix[MAXCELLS][j] = (float)0;		// Clearing the last column of the Matching_Matrix[][] for storing info on matched cells used below to assign new numbers to unmatched cells 

		// Repeatedly finding the lowest distance value in the  Matching_Matrix[][] and assuming it represents the same cell in this level and the previous. Continuing until the last matched cell movement is above the limit Max_Cell_Shift_px
		cells_matched = TRUE; 
		do{
			// Finding current best fit (=smallest value) in Matching_Matrix[][] 
			smallest_shift_1 = BIG;							// Just setting a large pixel distance 
			for( i=1; i<=Nof__Cells[level-1]; i++ )
				for( j=1; j<=Nof__Cells[level]; j++ )
					if( (shift=Matching_Matrix[i][j]) < smallest_shift_1 ) { smallest_shift_1=shift; i_1=i; j_1=j; } 

/*			// Finding the second nearest cell in this level to the potential mother cell in the previous level   
			smallest_shift_2 = BIG; 
			for( j=1; j<=Nof__Cells[level]; j++ )
					if( j!=j_1 && (shift=Matching_Matrix[i_1][j]) < smallest_shift_2 ) { smallest_shift_2=shift; j_2=j; } 

			// Looking first for a possible cell division event, (i.e. if the second nearest cell is similarly close, their areas are nearly half that of the mother cell, the two daughter cells lie in opposite directions from the mother cell )
			x_1 = Cells[level][j_1].cx-Cells[level-1][i_1].cx;		// Vector from mother cell to nearest cell in current level
			y_1 = Cells[level][j_1].cy-Cells[level-1][i_1].cy;
			x_2 = Cells[level][j_2].cx-Cells[level-1][i_1].cx;		// Vector from mother cell to second nearest cell in current level
			y_2 = Cells[level][j_2].cy-Cells[level-1][i_1].cy;
			if( smallest_shift_1>0 && smallest_shift_2>0 ) {					// Just to be able to check that the value is valid before making the acos() calculation since roundoff errors seem to occur 
				core = ( x_1*x_2 + y_1*y_2 )/( smallest_shift_1*smallest_shift_2 );	
				if( core >  1 ) core =  1; 
				if( core < -1 ) core = -1; 
				alpha = 180*acos(core)/Pi;										// Divergence angle [°] between vectors from mother cell to the daughter cell candidates 
			} else alpha = 999;													// Just a large value if cell has not moved so no shift vector to calculate angle from 

			// Special test of possible cell division NB: CURRENTLY NOT IN USE SO TEXT ABOVE IS COMMENTED OUT */
			if( FALSE && smallest_shift_1>0.3*Min_Diameter_px && smallest_shift_1<DIVISION_MARGIN*Min_Diameter_px && smallest_shift_2<DIVISION_MARGIN*Min_Diameter_px && alpha>DIVISION_ANGLE ){		

				for( i=1; i<=Nof__Cells[level-1]; i++ ) Matching_Matrix[i][j_1]=BIG;	Matching_Matrix[MAXCELLS][j_1]=BIG;		// Erasing the     best matched row in the Matching_Matrix[][] by setting to very high values to allow the nest best match to be found in the next iteration. Marking also the last column of the Matching_Matrix[][] for easy search below 
				for( i=1; i<=Nof__Cells[level-1]; i++ ) Matching_Matrix[i][j_2]=BIG;	Matching_Matrix[MAXCELLS][j_2]=BIG;		// Also the second best matched row in the Matching_Matrix[][] when one cell divides into two  
				for( j=1; j<=Nof__Cells[level  ]; j++ ) Matching_Matrix[i_1][j]=BIG;											// Erasing the column of the mother cell in the Matching_Matrix[][]  

				// Giving new cellnumber to the best fit daughter Cell (number j_1) and painting its pixels with the NEGATIVE of that number
				Nof_Matched_Cells++;															// Counting up the number of matched cells, since introducing a newborn cell  
				if( Nof_Matched_Cells == MATCHED_CELLS ) { sprintf(message, "The number of matched cells (Nof_Matched_Cells = %d) has reached\n the maximum size of the arrays which store cell data at image number %d.\nReduce the number of images or increase MATCHED_CELLS and recompile.\n", Nof_Matched_Cells, level+First_Image_Number);  write_error_message(ErrorMessageFile, ERROR, message); }
				Cells[level][j_1].cellnumber = Nof_Matched_Cells;								// Assigning the new cell with a new, unique number, which will then be assigned also to matching Cells[][] at other levels to allow drawing matched cells with the same color (since they will typically have different initial index value) 
				Matched_Cells[level][Nof_Matched_Cells] = j_1;									// Storing the original cell index at this level for the matched cell number in the array Matched_Cells[][] 
				fill_segmented(j_1, -Nof_Matched_Cells, Cells[level][j_1].cm_x_first, Cells[level][j_1].cm_y_first);	// Marking NEW cells in this level as "matched" for the next level by filling it with the NEGATIVE of its cell number 

				// Giving new cellnumber to the second best fit daughter Cell (number j_2) and painting its pixels with the NEGATIVE of that number
				Nof_Matched_Cells++;															// Counting up the number of matched cells, since introducing a newborn cell  
				if( Nof_Matched_Cells == MATCHED_CELLS ) { sprintf(message, "The number of matched cells (Nof_Matched_Cells = %d) has reached\n the maximum size of the arrays which store cell data at image number %d.\nReduce the number of images or increase MATCHED_CELLS and recompile.\n", Nof_Matched_Cells, level+First_Image_Number);  write_error_message(ErrorMessageFile, ERROR, message); }
				Cells[level][j_2].cellnumber = Nof_Matched_Cells;								// Assigning the new cell with a new, unique number, which will then be assigned also to matching Cells[][] at other levels to allow drawing matched cells with the same color (since they will typically have different initial index value) 
				Matched_Cells[level][Nof_Matched_Cells] = j_2;									// Storing the original cell index at this level for the matched cell number in the array Matched_Cells[][] 
				fill_segmented(j_2, -Nof_Matched_Cells, Cells[level][j_2].cm_x_first, Cells[level][j_2].cm_y_first);	// Marking NEW cells in this level as "matched" for the next level by filling it with the NEGATIVE of its cell number 

/*				// Temporarily drawing bounding boxes around all matched cells  
				for( x=Cells[level][j_1].min_x; x<=Cells[level][j_1].max_x; x++ ) Segmented_Image[SHIFTXY+x][SHIFTXY+Cells[level][j_1].min_y] = Segmented_Image[SHIFTXY+x][SHIFTXY+Cells[level][j_1].max_y] = -Cells[level][j_1].cellnumber;						
				for( y=Cells[level][j_1].min_y; y<=Cells[level][j_1].max_y; y++ ) Segmented_Image[SHIFTXY+Cells[level][j_1].min_x][SHIFTXY+y] = Segmented_Image[SHIFTXY+Cells[level][j_1].max_x][SHIFTXY+y] = -Cells[level][j_1].cellnumber;				
				for( x=Cells[level][j_2].min_x; x<=Cells[level][j_2].max_x; x++ ) Segmented_Image[SHIFTXY+x][SHIFTXY+Cells[level][j_2].min_y] = Segmented_Image[SHIFTXY+x][SHIFTXY+Cells[level][j_2].max_y] = -Cells[level][j_2].cellnumber;						
				for( y=Cells[level][j_2].min_y; y<=Cells[level][j_2].max_y; y++ ) Segmented_Image[SHIFTXY+Cells[level][j_2].min_x][SHIFTXY+y] = Segmented_Image[SHIFTXY+Cells[level][j_2].max_x][SHIFTXY+y] = -Cells[level][j_2].cellnumber;			//*/ 			

			// Looking for ordinary matching of one cell with itself in the next layer: Finding and matching the closest candidate until the given Max_Cell_Shift_px limitation is reached 
			} else if( smallest_shift_1 < Max_Cell_Shift_px ){	 

//				printf("\nlevel = %3d   i_1 = %4d   j_1 = %4d   smallest_shift_1 = %8.3f px = %8.3f um   Max_Cell_Shift_px = %8.3f um   Nof__Cells[level-1] = %d    Nof__Cells[level] = %d \n", level, i_1, j_1, smallest_shift_1, smallest_shift_1*Pixel_Size, Max_Cell_Shift_px*Pixel_Size, Nof__Cells[level-1], Nof__Cells[level]); 
				Matching_Matrix[MAXCELLS][j_1]=BIG;													// Marking the last column of the Matching_Matrix[][] for easy search below 
				for( i=1; i<=Nof__Cells[level-1]; i++ ) Matching_Matrix[i  ][j_1]=BIG;				// Erasing the matched rows and columns in the Matching_Matrix[][] by setting to very high values to allow the next best match to be found in the subsequent iteration. Marking also the column just right of the Matching_Matrix[][] for easy search below 
				for( j=1; j<=Nof__Cells[level  ]; j++ ) Matching_Matrix[i_1][  j]=BIG;

//				printf("Before Cells[level][j_1].cellnumber = Cells[level-1][i_1].cellnumber...");
				Cells[level][j_1].cellnumber = Cells[level-1][i_1].cellnumber;						// Propagating the unique cell identification number to the matched cell in the level above to allow drawing with the same color 
//				printf("Cells[%2d][%4d].cellnumber  = %4d   =   Cells[%2d][%4d].cellnumber  = %4d  \n", level, j_1, Cells[level][j_1].cellnumber, level-1, i_1, Cells[level-1][i_1].cellnumber);
				Matched_Cells[level][ Cells[level][j_1].cellnumber ] = j_1; 
				fill_segmented(j_1,  -Cells[level][j_1].cellnumber, Cells[level][j_1].cm_x_first, Cells[level][j_1].cm_y_first);		// Marking the matched cell in this level by filling them with the NEGATIVE of their cell number, which corresponds to the number of that cell in the previous level  

/*				// Temporarily drawing bounding boxes around all matched cells  
				for( x=Cells[level][j_1].min_x; x<=Cells[level][j_1].max_x; x++ ) Segmented_Image[SHIFTXY+x][SHIFTXY+Cells[level][j_1].min_y] = Segmented_Image[SHIFTXY+x][SHIFTXY+Cells[level][j_1].max_y] = -Cells[level][j_1].cellnumber;						
				for( y=Cells[level][j_1].min_y; y<=Cells[level][j_1].max_y; y++ ) Segmented_Image[SHIFTXY+Cells[level][j_1].min_x][SHIFTXY+y] = Segmented_Image[SHIFTXY+Cells[level][j_1].max_x][SHIFTXY+y] = -Cells[level][j_1].cellnumber;			//*/ 			

			} else cells_matched = FALSE; 

		}  while ( cells_matched );
//		printf("Finished do...while( cells_matched )\n"); 

		// For the remaining non-matched cells at this level, trying to connect track if cell is missing in previous image, otherwise assigning new numbers to the non-matched cells since they should be included in the matching with higher levels 
		for( j_1=1; j_1<=Nof__Cells[level]; j_1++ ){

			if( Matching_Matrix[MAXCELLS][j_1] < BIG-1 ){			// Info on which cells in this level that have been matches is stored in the last column of the Matching_Matrix[][]

				cells_matched = FALSE;

				// Searching TWO levels below for a possibly matching cell as a way of closing tracks that lack only one cell detection in the previous level 
				if( level > 1 ) { 

					smallest_shift_1 = BIG;							// Just setting a large pixel distance 
					for( i=1; i<=Nof__Cells[level-2]; i++ ){
						shift = (float)sqrt( (Cells[level][j_1].dx-Cells[level-2][i].dx)*(Cells[level][j_1].dx-Cells[level-2][i].dx) + (Cells[level][j_1].dy-Cells[level-2][i].dy)*(Cells[level][j_1].dy-Cells[level-2][i].dy) );
						if( shift < smallest_shift_1 ) { smallest_shift_1=shift; i_1=i; } 
					}

					if( smallest_shift_1 < 2*Max_Cell_Shift_px && Matched_Cells[level-1][ Cells[level-2][i_1].cellnumber ]==0 ){	// Found a cell two levels below with separation less twice the maximum step between two consecutive images and which is not matched to a cell one level below. 

						Cells[level][j_1].cellnumber = Cells[level-2][i_1].cellnumber;				// Propagating the unique cell identification number to the matched cell in the level two steps above to allow drawing with the same color 
						Matched_Cells[level][ Cells[level][j_1].cellnumber ] = j_1; 
						fill_segmented(j_1,  -Cells[level][j_1].cellnumber, Cells[level][j_1].cm_x_first, Cells[level][j_1].cm_y_first);		// Marking the matched cell in this level by filling them with the NEGATIVE of their cell number, which corresponds to the number of that cell in the previous level  

						// Since this cell was not detected in the previous level, adding a new cell to that level and interpolating cell properties as a best estimate
						j_2 = ++Nof__Cells[level-1];												// Reusing j_2 as shorthand 	
						if( j_2 == MAXCELLS-1 ) { sprintf(message, "The number of cells (Nof__Cells[level-1] = %d) has reached the \nsize of the arrays which store cell data (dimension = %d) at image number %d.\nReduce the number of images or increase MAXCELLS and recompile.\n", Nof__Cells[level-1], MAXCELLS, level+First_Image_Number);  write_error_message(ErrorMessageFile, ERROR, message); }

						Cells[level-1][j_2].cellnumber = Cells[level-2][i_1].cellnumber;			// Copying the unique cell identification number also to the new, interpolated cell in the intermediate level to allow drawing with the same color 
						Matched_Cells[level-1][ Cells[level-1][j_2].cellnumber ] = j_2; 

						Cells[level-1][j_2].area_fraction =       0.5*( Cells[level-2][i_1].area_fraction + Cells[level][j_1].area_fraction );	 
						Cells[level-1][j_2].aspect_ratio  =       0.5*( Cells[level-2][i_1].aspect_ratio  + Cells[level][j_1].aspect_ratio  );	 
						Cells[level-1][j_2].area_px       = (int)(0.5*( Cells[level-2][i_1].area_px       + Cells[level][j_1].area_px       ));	 
						Cells[level-1][j_2].cm_x_first    = (int)(0.5*( Cells[level-2][i_1].cm_x_first    + Cells[level][j_1].cm_x_first    ));	 
						Cells[level-1][j_2].cm_y_first    = (int)(0.5*( Cells[level-2][i_1].cm_y_first    + Cells[level][j_1].cm_y_first    ));	 
						Cells[level-1][j_2].cx            = (int)(0.5*( Cells[level-2][i_1].cx            + Cells[level][j_1].cx            ));	 
						Cells[level-1][j_2].cy            = (int)(0.5*( Cells[level-2][i_1].cy            + Cells[level][j_1].cy            ));	 
						Cells[level-1][j_2].dx            =       0.5*( Cells[level-2][i_1].dx            + Cells[level][j_1].dx            );	 
						Cells[level-1][j_2].dy            =       0.5*( Cells[level-2][i_1].dy            + Cells[level][j_1].dy            );	 
						Cells[level-1][j_2].fx            =       0.5*( Cells[level-2][i_1].fx            + Cells[level][j_1].fx            );	 
						Cells[level-1][j_2].fy            =       0.5*( Cells[level-2][i_1].fy            + Cells[level][j_1].fy            );	 
						Cells[level-1][j_2].max_x         = (int)(0.5*( Cells[level-2][i_1].max_x         + Cells[level][j_1].max_x         ));	 
						Cells[level-1][j_2].max_y         = (int)(0.5*( Cells[level-2][i_1].max_y         + Cells[level][j_1].max_y         ));	 
						Cells[level-1][j_2].min_x         = (int)(0.5*( Cells[level-2][i_1].min_x         + Cells[level][j_1].min_x         ));	 
						Cells[level-1][j_2].min_y         = (int)(0.5*( Cells[level-2][i_1].min_y         + Cells[level][j_1].min_y         ));		//*/
						cells_matched = TRUE; 
					}  
				}

				if( !cells_matched ){
					Nof_Matched_Cells++;														// Counting up the number of matched cells to introduce a new, non-matched cell  
					if( Nof_Matched_Cells == MATCHED_CELLS ) { sprintf(message, "The number of matched cells (Nof_Matched_Cells = %d) has reached\n the maximum size of the arrays which store cell data at image number %d.\nReduce the number of images or increase MATCHED_CELLS and recompile.\n", Nof_Matched_Cells, level+First_Image_Number);  write_error_message(ErrorMessageFile, ERROR, message); }
					Cells[level][j_1].cellnumber = Nof_Matched_Cells;							// Assigning each physical cell identified in the first image a unique number, which will then be assigned also to Cells[][] at other levels (which will typically have different last index value) 
					Matched_Cells[level][Nof_Matched_Cells] = j_1;								// Storing the original cell index at this level for the matched cell number in the array Matched_Cells[][] 
					fill_segmented(j_1, -Nof_Matched_Cells, Cells[level][j_1].cm_x_first, Cells[level][j_1].cm_y_first);	// Marking NEW cells in this level as "matched" for the next level by filling it with the NEGATIVE of its cell number 
				}
			}		
		}
		printf(          "%5d     %5d\n", Nof__Cells[level], Nof_Matched_Cells);
		fprintf(CMDWindowFile, "%d,%d\n", Nof__Cells[level], Nof_Matched_Cells);
		Nof_Identified_Cells[level] = Nof__Cells[level]; 
	} 

	if( Wound_Healing_Mode ) {

		// First zeroing, then populating histograms for number of matched cells at given y coordinate positions, first the distribution, then the accumulated distribution in both directions from y center of Bimage 
		for( y=0; y<Bimage->ysize; y++ ) WH_Matched_Histogram[0][y] = WH_Matched_Histogram[1][y] = 0; 
		for( j=1; j<=Nof_Matched_Cells; j++ ){
			y = Cells[level][Matched_Cells[level][j]].cy;
			if( y>0 && y<MAXY ) WH_Matched_Histogram[0][y]++;
		}

		y_1 = (int)(0.5*Bimage->ysize);													// Just reusing variable as midpoint of cropped image, assumed to be the center of the wound 
		if( level == 0 ) WH_Nof_Cells_Upper_First = WH_Nof_Cells_Lower_First = BIG;		// Setting large limits for the zeroth image to find the proper starting limits 

		// Counting number of matched cells in LOWER half of image and finding the position (WH_Box_Lower_Limit_px) where the number of cells matches the cellcount in the lower half of the first image 
		for( y=y_1-1; y>0; y-- ) WH_Matched_Histogram[1][y] = WH_Matched_Histogram[1][y+1] + WH_Matched_Histogram[0][y];
		if( WH_Matched_Histogram[1][1] > WH_Nof_Cells_Lower_First ) {
			y = 0;
			while( WH_Matched_Histogram[1][++y] > WH_Nof_Cells_Lower_First );
			WH_Box_Lower_Limit_px = y;													// Finding the lower edge in the cropped image with the original number of cells above it 
		} else WH_Box_Lower_Limit_px -= 3;												// Assuming 3 pixels back movement is enough in most videos. Looked slightly better than 1 in a test.  
		if( WH_Box_Lower_Limit_px < 0 ) WH_Box_Lower_Limit_px = 0;  

		// Counting number of matched cells in UPPER half of image and finding the position (WH_Box_Upper_Limit_px) where the number of cells matches the cellcount in the lower half of the first image 
		for( y=y_1+1; y<Bimage->ysize; y++ ) WH_Matched_Histogram[1][y] = WH_Matched_Histogram[1][y-1] + WH_Matched_Histogram[0][y];
		if( WH_Matched_Histogram[1][Bimage->ysize-1] > WH_Nof_Cells_Upper_First ) {
			y = Bimage->ysize;
			while( WH_Matched_Histogram[1][--y] > WH_Nof_Cells_Upper_First );
			WH_Box_Upper_Limit_px = y;													// Finding the upper edge in the cropped image with the original number of cells below it 
		} else WH_Box_Upper_Limit_px += 3;
		if( WH_Box_Upper_Limit_px >= Bimage->ysize ) WH_Box_Upper_Limit_px = Bimage->ysize-1;  

	/*	for( y=0; y<=Bimage->ysize; y++ ) {
			fprintf(CMDWindowFile, "%d,%d,%d,%d,%d,%d,%d\n", y, WH_Matched_Histogram[0][y], WH_Matched_Histogram[1][y], WH_Box_Lower_Limit_px, WH_Box_Upper_Limit_px, WH_Nof_Cells_Lower, WH_Nof_Cells_Upper);
		} //*/
		if( level == 0) {													// After analyzing the zeroth image, setting the cell counts to the actual values  
			WH_Nof_Cells_Lower_First = WH_Matched_Histogram[1][1];
			WH_Nof_Cells_Upper_First = WH_Matched_Histogram[1][Bimage->ysize-1];
			// printf("Box limit: %6d %6d   Nof cells %6d %6d  \n", WH_Box_Lower_Limit_px, WH_Box_Upper_Limit_px, WH_Nof_Cells_Lower_First, WH_Nof_Cells_Upper_First);
		}
	}

	// Writing MATCHED image to file for checking 
	copy_bitmap(Rimage, Oimage);											// Copying original image into output image 
	convert2color(Oimage);													// To allow showing different cells with different colors 

	// Drawing the track of all cells that are currently matched as a line of the same color as the cell, to see more easily how much the cell has moved 
	if( Draw_Track && level > 0 ){
		for( j=1; j<=Nof_Matched_Cells; j++ ){
			m=Cells[level][Matched_Cells[level][j]].cellnumber;							// The cell number used to get the same color as the cell in the image 
			for( lev=level; lev>0; lev-- ){
				if( Matched_Cells[level][j]>0 && (m_2=Matched_Cells[lev][j])>0 && (m_1=Matched_Cells[lev-1][j])>0 ){					// We have a matched cell in this level and also the two current levels for drawing  
					x_2 = Cells[lev  ][m_2].cx;											// Shorthand 
					y_2 = Cells[lev  ][m_2].cy;											// Unit: [pixels]
					x_1 = Cells[lev-1][m_1].cx;			
					y_1 = Cells[lev-1][m_1].cy;
					delta_x = x_2 - x_1; 
					delta_y = y_2 - y_1; 
					delta_max = maxi(abs(delta_x), abs(delta_y));
					if( delta_max == 0 ) delta_max = 1; 
					for( d=0; d<=delta_max; d++ ) {
						cellnumber2color(m, color);
						for( i=0; i<3; i++ ) Oimage->image[(int)(Left___Crop_Margin+x_1+1.0*d*delta_x/delta_max+0.5)][(int)(Bottom_Crop_Margin+y_1+1.0*d*delta_y/delta_max+0.5)][i] = color[i]; 	
					}
				}
			}
		}
	} 

	// Drawing a filled circle with the correct cell color at the center of each cell in the original image 
	dot_radius = (float)(0.5*Dot_Diameter_um/Pixel_Size);		// Radius in um  
	int_radius = (int)(dot_radius+0.5);							// Rounding off the diameter to find integer limits to the two loops below  
	if( int_radius < 1 ) int_radius = 1;						// Using at least 3 (=1+2*int_radius) as smallest spot diameter to draw in original image at cell center of mass
	if( Draw_CM ) 
		for( j=1; j<=Nof__Cells[level]; j++ )
			for( y=-int_radius; y<=int_radius; y++ )
				for( x=-int_radius; x<=int_radius; x++ )
					if( sqrt(1.0*x*x+y*y)<dot_radius && Left___Crop_Margin+Cells[level][j].cx+x>0 && Bottom_Crop_Margin+Cells[level][j].cy+y>0 && Left___Crop_Margin+Cells[level][j].cx+x<Rimage->xsize && Bottom_Crop_Margin+Cells[level][j].cy+y<Rimage->ysize ) {
						cellnumber2color(Cells[level][j].cellnumber,color); 
						for( i=0; i<3; i++ ) Oimage->image[Left___Crop_Margin+Cells[level][j].cx+x][Bottom_Crop_Margin+Cells[level][j].cy+y][i] = color[i];
					}

	// Drawing the contours of matched cells onto original image 
	if( Draw_Outline ) 
		for( y=0; y<Bimage->ysize; y++ )
			for( x=0; x<Bimage->xsize; x++ ) 
				if( (s=Segmented_Image[SHIFTXY+x][SHIFTXY+y])<0 && (Segmented_Image[SHIFTXY+x+1][SHIFTXY+y]!=s || Segmented_Image[SHIFTXY+x-1][SHIFTXY+y]!=s || Segmented_Image[SHIFTXY+x][SHIFTXY+y+1]!=s || Segmented_Image[SHIFTXY+x][SHIFTXY+y-1]!=s ) ) {
					cellnumber2color(-s,color);
					for( i=0; i<3; i++ ) Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][i] = color[i]; 			//*/
				}

	// Drawing horizontal, white lines corresponding to the lower and upper edge of the current wound healing inclusion box  
	if( Wound_Healing_Mode ) {
		for( x=0; x<Bimage->xsize; x+=2 ) {
			Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+WH_Box_Lower_Limit_px][0] = Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+WH_Box_Lower_Limit_px][1] = Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+WH_Box_Lower_Limit_px][2] = (unsigned char)255;
			Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+WH_Box_Upper_Limit_px][0] = Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+WH_Box_Upper_Limit_px][1] = Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+WH_Box_Upper_Limit_px][2] = (unsigned char)255;
		}
	}

/*	// Drawing a single YELLOW pixel at the center of mass of cells in the previous level, to see more easily how much cells have moved 
	if( level > 0 )
		for( j=1; j<=Nof__Cells[level-1]; j++ ){
			Oimage->image[Left___Crop_Margin+Cells[level-1][j].cx][Bottom_Crop_Margin+Cells[level-1][j].cy][0] =   0; 	
			Oimage->image[Left___Crop_Margin+Cells[level-1][j].cx][Bottom_Crop_Margin+Cells[level-1][j].cy][1] = 255; 	
			Oimage->image[Left___Crop_Margin+Cells[level-1][j].cx][Bottom_Crop_Margin+Cells[level-1][j].cy][2] = 255; 	
		}		//*/

/*	// Drawing the treated image into the original 
	for( y=0; y<Bimage->ysize; y++ )
		for( x=0; x<Bimage->xsize; x++ ) {
			cellnumber2color(-Segmented_Image[SHIFTXY+x][SHIFTXY+y],color); 	
			for( i=0; i<3; i++) Oimage->image[y+Left___Crop_Margin][y+Bottom_Crop_Margin][i] = color[i]; 	
		}	//*/

	// Writing MATCHED image to file for making video or for checking 
	if( Write_6_MA_img || Make__6_MA_vid ) {
		if( Add_Scalebar ) draw_scalebar(); 
		sprintf(filename, "%s/%s_06_%04d_Matched.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*level);
		out_bitmap(Oimage, filename);
	}
}


// Smoothing the cell tracks over time by a simple gliding average filtering with weights 1:2:1 for current time step (level) and the levels before and after. Keeping the first and last track position fixed.

void smooth_tracks(int level)
{
	int i, n, lev, p_nr, t_nr, n_nr;
	double cm_x_prev, cm_y_prev, cm_x_this, cm_y_this, cm_x_next, cm_y_next; 

	printf("Smoothing cell track, using %d iterations.\n", Track_Smoothing_Iterations); 

	for( n=0; n<Track_Smoothing_Iterations; n++ ){ 

		// Copying all fractional cell positions [dx,dy] to temporary variables [fx,fy] in preparation for cell track filtering 
		for( i=1; i<Nof_Matched_Cells; i++ ){
			for( lev=0; lev<level; lev++ ){
				Cells[lev][i].fx = Cells[lev][i].dx;
				Cells[lev][i].fy = Cells[lev][i].dy;
			}
		}

		for( i=1; i<Nof_Matched_Cells; i++ ){
			for( lev=1; lev<level-1; lev++ ){

				if( (p_nr=Matched_Cells[lev-1][i]) > 0 ){										
					cm_x_prev = Cells[lev-1][p_nr].fx;											// Shorthand + giving the possibility to mark as -1 if not a tack 
					cm_y_prev = Cells[lev-1][p_nr].fy;				
				} else cm_x_prev = cm_y_prev = -1;				

				if( (t_nr=Matched_Cells[lev  ][i]) > 0 ){										// Calculating center of mass coordinates only for valid matches 
					cm_x_this = Cells[lev  ][t_nr].fx;											 
					cm_y_this = Cells[lev  ][t_nr].fy;											// Unit: [pixels]
				} else cm_x_this = cm_y_this = -1;												// Erasing old value if no match 

				if( (n_nr=Matched_Cells[lev+1][i]) > 0 ){
					cm_x_next = Cells[lev+1][n_nr].fx;
					cm_y_next = Cells[lev+1][n_nr].fy;
				} else cm_x_next = cm_y_next = -1;		

				if( cm_x_prev>0 && cm_x_this>0 && cm_x_next>0 ){								// If three consecutive valid values, calculating special gliding average with 1:2:1 weights which will flatten an alternating signal 
					Cells[lev][t_nr].dx = 0.25*cm_x_prev + 0.5*cm_x_this + 0.25*cm_x_next;		// Writing back to variables dx and dy which are used for calculating cell velocities etc. later		
					Cells[lev][t_nr].dy = 0.25*cm_y_prev + 0.5*cm_y_this + 0.25*cm_y_next;	
					Cells[lev][t_nr].cx = (int)( Cells[lev][t_nr].dx+0.5 );						// Also updating integer coordinates for center-of-mass to be used when drawing valid tracks images below 		
					Cells[lev][t_nr].cy = (int)( Cells[lev][t_nr].dy+0.5 );								
				}
			}
		}
	}
}


// Calculates statistical results based on valid tracks and writes resykts to various csv-files

void finish_statistics_files(int level)
{
	char ch, well_letter, video_number_string[9];
	unsigned char color[3]; 
	unsigned int l;
	int i, j, lev, c_nr, hist_index, velocity_histmax, video_number, well_number;
	int first_level, last__level, nof_trax_used;
	long m;
	double cm_x_this, cm_y_this, cm_x_next, cm_y_next; 
	double cm_x_first, cm_y_first, cm_x__last, cm_y__last; 
	double delta_x, delta_y, distance, angle; 
	double histogram_velocity_mean, histogram_velocity_mode__maximum, histogram_velocity_mode_parabola, histogram_velocity_median, histogram_velocity_75_percentile; 
	double             velocity,             velocity_sum,             velocity_sum2; 
	double accumulated_distance, accumulated_distance_sum, accumulated_distance_sum2, accumulated_distance_mean, accumulated_distance_sd, accumulated_distance_sem;  
	double   euclidean_distance,   euclidean_distance_sum,   euclidean_distance_sum2,   euclidean_distance_mean,   euclidean_distance_sd,   euclidean_distance_sem; 
	double   euclidean_velocity,   euclidean_velocity_sum,   euclidean_velocity_sum2,   euclidean_velocity_mean,   euclidean_velocity_sd,   euclidean_velocity_sem;
	double           directness,           directness_sum,           direcntess_sum2,           directness_mean,           directness_sd,           directness_sem;
	double                fmi_x,                fmi_x_sum,                fmi_x_sum2,                fmi_x_mean,                fmi_x_sd,                fmi_x_sem; 
	double                fmi_y,                fmi_y_sum,                fmi_y_sum2,                fmi_y_mean,                fmi_y_sd,                fmi_y_sem;
	// Parabola fitting parameters 											
	int fit_start, fit_midpoint;
	int n, halflength;									 
	double a, c, e, p;							// Coefficients for least square fit parabola to 2p+1 points in profile 
	double z0, z1, z2;							// Sum values for right hand side in matrix equation                    
	double A, B, C;								// Coefficients in z = Axx + Bx + C decribing the fitting parabolas		 

	// Going through the Matched_Cells[][] matrix and counting the number of images over which each cell has been matched 
	for( i=1; i<=Nof_Matched_Cells; i++ ){
		Matched_Cells[level][i]=0;						// Using the bottom row in the Matched_Cells[][] matrix to store the number of levels over which this cell has been matched. Done to allow outputting data only for cells matched more than Shortest_Cell_Track times. This removes dubious cells and noise. 
		for( lev=0; lev<level; lev++ )
			if( Matched_Cells[lev  ][i] > 0 )
				Matched_Cells[level][i]++;
	}

	// Writing the full Matched_Cells[][] matrix to file for examination 
	fprintf(MatrixFile, "Level / Match #,");
	for( i=1; i<=Nof_Matched_Cells; i++ ) fprintf(MatrixFile, "%d,", i);
	fprintf(MatrixFile, "\n");
	for( lev=0; lev<level; lev++ ){
		fprintf(MatrixFile, "%d,", lev);
		for( i=1; i<=Nof_Matched_Cells; i++ ) fprintf(MatrixFile, "%d,", Matched_Cells[lev][i]);
		fprintf(MatrixFile, "\n");
	} 
	fprintf(MatrixFile, "# matched,");
	for( i=1; i<=Nof_Matched_Cells; i++ ) fprintf(MatrixFile, "%d,", Matched_Cells[level][i]);
	fprintf(MatrixFile, "\n");			//*/

	// Writing position and track data to files, and validated images containing only valid tracks
	fprintf(PosFile, "Track number,Image number,Time [min],Image shift X [um],Image shift Y [um],Cell center X [um],Cell center Y [um],Cell shift X [um],Cell shift Y [um],Distance [um],Velocity [um/min],Area [pixels],Aspect ratio [1],Box fill [1],Hex RGB color\n");	
//  POSSIBLY ADD LATER:        START TIME [h],END TIME [h]         
	fprintf(TraxFile,"Track number,Start image,End image,Start X [um],Start Y [um],End X [um],End Y [um],Delta X [um],Delta Y [um],Euclidean distance [um],Euclidean velocity [um/min],Accumulated distance [um],Directness [Eucl/Accum],Forward migration index X [Delta X/Accum],Forward migration index Y [Delta Y/Accum],Angle [°],Number of time steps,Mean velocity [um/min],SD velocity [um/min],SEM velocity [um/min]\n"); 
	nof_trax_used=0;			// Counter for matched cells found over more than Shortest_Cell_Track levels, to be used in Velocity_Matrix[][]  
	accumulated_distance_sum = accumulated_distance_sum2 = euclidean_distance_sum = euclidean_distance_sum2 = euclidean_velocity_sum = euclidean_velocity_sum2 = directness_sum = direcntess_sum2 = fmi_x_sum = fmi_x_sum2 = fmi_y_sum = fmi_y_sum2 = 0;  // Just to be sure

	for( i=1; i<=Nof_Matched_Cells; i++ ){	
		if( Matched_Cells[level][i] >= Shortest_Cell_Track ){			// Only processing tracks that have been matched for at least Shortest_Cell_Track images	
			nof_trax_used++;											// Counting up number of cells to store in Velocity_Matrix[][]
			first_level = last__level = -1;								// Clearing the variables for the first and last level a cell has been matched in 
			accumulated_distance = velocity_sum = velocity_sum2 = 0;	// Zeroing variables to be summed over all cell movements below 

//			for( lev=0; lev<level; lev++ ){
			for( lev=0; lev<=level; lev++ ){

				cellnumber2color(Cells[lev][Matched_Cells[lev][i]].cellnumber, color);						// To allow writing the track color in the Position data file

				if( (c_nr=Matched_Cells[lev  ][i]) > 0 ){													// Calculating center of mass coordinates only for valid matches 
					cm_x_this = Cells[lev  ][c_nr].dx*Pixel_Size;											// Shorthand 
					cm_y_this = Cells[lev  ][c_nr].dy*Pixel_Size;											// Unit: [um]
				} else cm_x_this = cm_y_this = -1;															// Erasing old value if no match 

				if( (c_nr=Matched_Cells[lev+1][i]) > 0 ){
					cm_x_next = Cells[lev+1][c_nr].dx*Pixel_Size;
					cm_y_next = Cells[lev+1][c_nr].dy*Pixel_Size;
				} else cm_x_next = cm_y_next = -1;

				if( first_level==-1 && Matched_Cells[lev][i]>0 ){
					first_level = lev;																		// The first level with a match in the Matched_Cells[][] matrix 
					cm_x_first = cm_x_this;																	// Shorthand 
					cm_y_first = cm_y_this;																	// Unit: [um]		
				}
				if( last__level==-1 && (lev==level-1 || (first_level>=0 && Matched_Cells[lev+1][i]==0 ))){
					last__level = lev;																		// The last level with a match in the Matched_Cells[][] matrix 
					cm_x__last = cm_x_this;			
					cm_y__last = cm_y_this;
				}

				if( cm_x_this>0 && cm_y_this>0 && cm_x_next>0 && cm_y_next>0 ) {
					distance = sqrt( (cm_x_next-cm_x_this)*(cm_x_next-cm_x_this)+(cm_y_next-cm_y_this)*(cm_y_next-cm_y_this) );			// Unit: [µm]
					accumulated_distance += distance; 													// Unit: [µm]
					velocity = distance/(TimeStep*Image_Nr_Increment);									// Unit: [µm/min]
					velocity_sum  += velocity;															// For mean and SD calculations
					velocity_sum2 += velocity*velocity;													// For mean and SD calculations
					hist_index = (int)( Nof_Bins*velocity/Max_Cell_Velocity );
					if( hist_index > Nof_Bins ) hist_index = Nof_Bins;  
					Velocity_Histogram[lev  ][hist_index]++;											// The level specific part of the Velocity_Histogram[][]
					Velocity_Histogram[level][hist_index]++; 											// The complete Velocity_Histogram[][] accumulated over all levels is stored here 
					Velocity_Matrix[lev][5+nof_trax_used] = (float)velocity;							// Putting velocity value at correct location in Velocity_Matrix[][]
					Cell_Shift_Vector_Sum_X[lev] += (float)(cm_x_next-cm_x_this);						// Adding X coordinate of total cell shift to calculate accumulated, average velocity vector below 
					Cell_Shift_Vector_Sum_Y[lev] += (float)(cm_y_next-cm_y_this);						// Adding Y coordinate of total cell shift to calculate accumulated, average velocity vector below 
					if( Wound_Healing_Mode ) {
						if( cm_y_this/Pixel_Size > 0.5*Bimage->ysize ) {								// Using pixel coordinates (not µm) 
							Cell_Shift_Vector_Sum_X_WH_Upper[lev] += (float)(cm_x_next-cm_x_this);		// Accumulating cell movement coordinates separately for cells in upper half of wound healing image
							Cell_Shift_Vector_Sum_Y_WH_Upper[lev] += (float)(cm_y_next-cm_y_this);	
							WH_Nof_Cells_Upper[lev]++;													// Counting exact number of cell movements to get correct average values in velocity matrix
						} else {
							Cell_Shift_Vector_Sum_X_WH_Lower[lev] += (float)(cm_x_next-cm_x_this);		// Accumulating cell movement coordinates separately for cells in lower half of wound healing image
							Cell_Shift_Vector_Sum_Y_WH_Lower[lev] += (float)(cm_y_next-cm_y_this);	
							WH_Nof_Cells_Lower[lev]++;
						}
					}																																																											// New 2023-09-10: Adding crop margin distance from lower left corner to cell center of mass to make cell coordinates independend of the crop margins and given relative to lower left corner of image 
					fprintf(PosFile, "%d,%d,%1.0f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%7.5f,%d,%7.5f,%7.5f,%06X\n", nof_trax_used, First_Image_Number+Image_Nr_Increment*lev, TimeStep*Image_Nr_Increment*lev, Image_Shift_um_X[lev], Image_Shift_um_Y[lev], cm_x_this+Left___Crop_Margin*Pixel_Size, cm_y_this+Bottom_Crop_Margin*Pixel_Size, cm_x_this-cm_x_first, cm_y_this-cm_y_first, distance, velocity, Cells[lev][Matched_Cells[lev][i]].area_px, Cells[lev][Matched_Cells[lev][i]].aspect_ratio, Cells[lev][Matched_Cells[lev][i]].area_fraction, 65536*color[2]+256*color[1]+color[0] );
//					fprintf(PosFile, "%d,%d,%1.0f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%7.5f,%d,%7.5f,%7.5f,%06X\n", nof_trax_used, First_Image_Number+Image_Nr_Increment*lev, TimeStep*Image_Nr_Increment*lev, Image_Shift_um_X[lev], Image_Shift_um_Y[lev], cm_x_this, cm_y_this, cm_x_this-cm_x_first, cm_y_this-cm_y_first, distance, velocity, Cells[lev][Matched_Cells[lev][i]].area_px, Cells[lev][Matched_Cells[lev][i]].aspect_ratio, Cells[lev][Matched_Cells[lev][i]].area_fraction, 65536*color[2]+256*color[1]+color[0] );
//				} else if( cm_x_this>0 && cm_y_this>0 && cm_x_next<0 && cm_y_next<0 ) {		// The last image analysed which has no cell shift, distance or velocity to the next image, but has cell center coordinates etc. is treated separately 
				} else if( cm_x_this>0 && cm_y_this>0  ) {		// The last image analysed which has no distance or velocity to the next image, but has cell center coordinates etc. is treated separately 
					fprintf(PosFile, "%d,%d,%1.0f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,,,%d,%7.5f,%7.5f,%06X\n", nof_trax_used, First_Image_Number+Image_Nr_Increment*lev, TimeStep*Image_Nr_Increment*lev, Image_Shift_um_X[lev], Image_Shift_um_Y[lev], cm_x_this+Left___Crop_Margin*Pixel_Size, cm_y_this+Bottom_Crop_Margin*Pixel_Size,  cm_x_this-cm_x_first, cm_y_this-cm_y_first, Cells[lev][Matched_Cells[lev][i]].area_px, Cells[lev][Matched_Cells[lev][i]].aspect_ratio, Cells[lev][Matched_Cells[lev][i]].area_fraction, 65536*color[2]+256*color[1]+color[0] );
					   Velocity_Matrix[lev][5+nof_trax_used] = -1;							// Marking non valid velocities as negative in Velocity_Matrix[][] to allow empty cells (not zero) in Excel
				} else Velocity_Matrix[lev][5+nof_trax_used] = -1;							// Marking non valid velocities as negative in Velocity_Matrix[][] to allow empty cells (not zero) in Excel

			}

			delta_x = cm_x__last-cm_x_first;					// Shorthand							
			delta_y = cm_y__last-cm_y_first;					// Unit: [um]						

			if( delta_x==0 ) angle = -999;		// Unit: [degrees] 
			else			 angle = atan2(delta_y, delta_x)*180/Pi;
			
			// Preparing calculations over all tracks of Mean, SD and SEM for these variables reported in the _Overall_summary.csv file 
			accumulated_distance_sum  += accumulated_distance;															// The accumulated_distance calculated in the above loop, is directly proportional to velocity_sum
			accumulated_distance_sum2 += accumulated_distance*accumulated_distance;

			euclidean_distance       = sqrt( delta_x*delta_x + delta_y*delta_y  );										// The straight line from start to end position for this cell
			euclidean_distance_sum  += euclidean_distance;
			euclidean_distance_sum2 += euclidean_distance*euclidean_distance;

			euclidean_velocity       = euclidean_distance/(TimeStep*Image_Nr_Increment*(last__level-first_level));		// Average velocity over which the cell was tracked
			euclidean_velocity_sum  += euclidean_velocity;
			euclidean_velocity_sum2 += euclidean_velocity*euclidean_velocity;

			directness       = euclidean_distance/accumulated_distance;													// Will be 1 if cell moves in a straight line
			directness_sum  += directness;
			direcntess_sum2 += directness*directness;

			fmi_x       = absd(delta_x)/accumulated_distance;															// Taking the absolute value to make the mean FMI always positive, useful for wound healing where opposing directions would otherwise cancel each other
			fmi_x_sum  += fmi_x;
			fmi_x_sum2 += fmi_x*fmi_x;
			
			fmi_y       = absd(delta_y)/accumulated_distance;
			fmi_y_sum  += fmi_y;
			fmi_y_sum2 += fmi_y*fmi_y;

			m = last__level-first_level;

			// Summary info stored below the Velocity_Matrix[][] and shifted 5 steps to the right:
			Velocity_Matrix[level  ][5+nof_trax_used] = (float)(first_level+First_Image_Number);								// First level where cell was matched 
			Velocity_Matrix[level+1][5+nof_trax_used] = (float)(last__level+First_Image_Number);								// Last  level where cell was matched
			Velocity_Matrix[level+2][5+nof_trax_used] = (float)(last__level-first_level);										// Number of levels in this track 
			Velocity_Matrix[level+3][5+nof_trax_used] = (float)(velocity_sum/m);												// Mean velocity for this track 
			Velocity_Matrix[level+4][5+nof_trax_used] = (float)sqrt( (velocity_sum2-velocity_sum*velocity_sum/m)/(m-1) );		// SD  of the velocity for this track 
			Velocity_Matrix[level+5][5+nof_trax_used] = (float)(Velocity_Matrix[level+4][5+nof_trax_used]/sqrt((double)m));		// SEM of the velocity for this track 
			Velocity_Matrix[level+6][5+nof_trax_used] = (float)(euclidean_velocity);											// Less affected by noise from jumping centers of stationary cells  

			fprintf(PosFile, "\n");			// One blank line after each cell to get abrupted line in Excel plot
			fprintf(TraxFile,"%d,%d,%d,%3.1f,%3.1f,%3.1f,%3.1f,%3.1f,%3.1f,%3.1f,%5.3f,%3.1f,%7.5f,%7.5f,%7.5f,%3.1f,%d,%7.5f,%7.5f,%7.5f\n", nof_trax_used, First_Image_Number+Image_Nr_Increment*first_level, First_Image_Number+Image_Nr_Increment*last__level, cm_x_first, cm_y_first, cm_x__last, cm_y__last, delta_x, delta_y, euclidean_distance, euclidean_velocity, accumulated_distance, directness, fmi_x, fmi_y, angle, m, Velocity_Matrix[level+3][5+nof_trax_used], Velocity_Matrix[level+4][5+nof_trax_used], Velocity_Matrix[level+5][5+nof_trax_used] ); 
		}
	} // End of for( i over Nof_Matched_Cells...) 

	// Finding the highest point in the overall curve Velocity_Histogram[level][], fitting a parabola around this point and finding histogram_velocity_mode__maximum as the top point 
	velocity_histmax=0;
	for( i=0; i<=Nof_Bins; i++ ){
		if( Velocity_Histogram[level][i] > velocity_histmax ) { 
			velocity_histmax = Velocity_Histogram[level][i]; 
			fit_midpoint=i;
		}	
	}
	fit_start = 1;											// Beginning of fit region set to second datapoint in histogram since the first tends to be abnormally high
	halflength = fit_midpoint-fit_start;					// Distance on each side from fit_midpoint to fit parabola
	p = (double)halflength;									// Shorthand since used a lot below 
	a = 2*(p*p*p*p*p/5 + p*p*p*p/2 + p*p*p/3 - p/30);		// Calculating the matrix coefficients 
	c = p*(p+1)*(2*p+1)/3;
	e = 2*p+1;
	z0=0; for( n=-halflength; n<=halflength; n++ ) z0 +=     Velocity_Histogram[level][fit_midpoint+n];	
	z1=0; for( n=-halflength; n<=halflength; n++ ) z1 +=   n*Velocity_Histogram[level][fit_midpoint+n]; 
	z2=0; for( n=-halflength; n<=halflength; n++ ) z2 += n*n*Velocity_Histogram[level][fit_midpoint+n]; 

	// Calculating coefficients for fitted parabola and finding the parabola top point if the fit is successful 
	if( c*c-a*e!=0 && (A=( c*z0 - e*z2 )/( c*c - a*e ))< 0 ) {																// A non-infinite and negative quadratic coefficient is reqiured for a valid fit 
		B = z1/c;
		C = ( c*z2 - a*z0 )/( c*c - a*e );		
		for( n=-halflength; n<=halflength; n++ ) Velocity_Histogram[level+2][fit_midpoint+n] = (int)(A*n*n + B*n + C);		// Calculating shape and storing in Velocity_Histogram[level+2][] for later output
		histogram_velocity_mode_parabola = (fit_midpoint-0.5*B/A+0.5)*Max_Cell_Velocity/Nof_Bins;							// Parabola has maximum value at -B/2A, but the origin is at fit_midpoint. Adding 0.5 to match the bin midpoint as used below. This position is used as an estimate for the mode value which is less sensitive to noise in the histogram
	} else {
		for( n=-halflength; n<=halflength; n++ ) Velocity_Histogram[level+2][fit_midpoint+n] = (int)(-1);					// Just marking curve as constant negative to show it is not valid 
		histogram_velocity_mode_parabola = (fit_midpoint+0.5)*Max_Cell_Velocity/Nof_Bins;									// Taking the best estimate as the midpoint of the highest histogram column  
	}

	// Writing velocity HISTOGRAM to file, while calculating the mean, mode and median velocity for all cell movements 
	fprintf(VelHistFile,"Bin midpoint [um/min],Total counts,Fitted parabola,");
	for( lev=0; lev<level; lev++ ) fprintf(VelHistFile,"%1.0f min,", lev*TimeStep*Image_Nr_Increment); 
	fprintf(VelHistFile,"\n"); 
	velocity_sum=velocity_histmax=0;
	for( i=0; i<=Nof_Bins; i++ ){
		fprintf(VelHistFile,"%7.5f,%d,", (i+0.5)*Max_Cell_Velocity/Nof_Bins, Velocity_Histogram[level][i]);
		if( Velocity_Histogram[level+2][i]>0 ) fprintf(VelHistFile,"%d,", Velocity_Histogram[level+2][i]); else fprintf(VelHistFile,",");
		for( lev=0; lev<level; lev++ ) fprintf(VelHistFile,"%d,", Velocity_Histogram[lev][i]); 
		fprintf(VelHistFile,"\n"); 
		if( Velocity_Histogram[level][i] > velocity_histmax ) { velocity_histmax = Velocity_Histogram[level][i]; histogram_velocity_mode__maximum = (i+0.5)*Max_Cell_Velocity/Nof_Bins; }		

		if( i == 0) Velocity_Histogram[level+1][i] = Velocity_Histogram[level][i];
		else		Velocity_Histogram[level+1][i] = Velocity_Histogram[level][i] + Velocity_Histogram[level+1][i-1];		// Storing the accumulated number of counts in Velocity_Histogram[level+1][] to facilitate finding the median value easily below
		velocity_sum += Velocity_Histogram[level][i]*(i+0.5)*Max_Cell_Velocity/Nof_Bins;
	}

	if( Velocity_Histogram[level+1][Nof_Bins] > 0 ) {
		histogram_velocity_mean = velocity_sum/Velocity_Histogram[level+1][Nof_Bins];										// Velocity_Histogram[level+1][Nof_Bins] = total number of counts 
		i=0;
		while( i<Nof_Bins && Velocity_Histogram[level+1][i]<0.5*Velocity_Histogram[level+1][Nof_Bins] ) i++;				// Searching the accumulated number of counts until more than half of the counts reached 
		if( Velocity_Histogram[level+1][i] == Velocity_Histogram[level+1][i-1] ) histogram_velocity_median = (Max_Cell_Velocity/Nof_Bins)*(i - 0.5);			// Should not happen unless double peaked distribution with zero section in center 
		else 																	 histogram_velocity_median = (Max_Cell_Velocity/Nof_Bins)*(i - 0.5 + (0.5*Velocity_Histogram[level+1][Nof_Bins]-Velocity_Histogram[level+1][i-1])/(1.0*Velocity_Histogram[level+1][i]-Velocity_Histogram[level+1][i-1]) );			// Interpolating to find median value with sub-bin precision 
		while( i<Nof_Bins && Velocity_Histogram[level+1][i]<0.75*Velocity_Histogram[level+1][Nof_Bins] ) i++;				// Searching the accumulated number of counts until more than 75 % of the counts reached 
		if( Velocity_Histogram[level+1][i] == Velocity_Histogram[level+1][i-1] ) histogram_velocity_75_percentile = (Max_Cell_Velocity/Nof_Bins)*(i - 0.5);			// Should not happen unless double peaked distribution with zero section in center 
		else 																	 histogram_velocity_75_percentile = (Max_Cell_Velocity/Nof_Bins)*(i - 0.5 + (0.75*Velocity_Histogram[level+1][Nof_Bins]-Velocity_Histogram[level+1][i-1])/(1.0*Velocity_Histogram[level+1][i]-Velocity_Histogram[level+1][i-1]) );			// Interpolating to find median value with sub-bin precision 
	} else histogram_velocity_mean = histogram_velocity_median = -1;														// Marking as negative if no movements were counted. (Should not happen.) 

	fprintf(VelHistFile, "Velocity mode (highest point),,%7.5f,um/min\n", histogram_velocity_mode__maximum );
	fprintf(VelHistFile, "Velocity mode (parabola fit),,%7.5f,um/min\n",  histogram_velocity_mode_parabola );
	fprintf(VelHistFile, "Velocity median,,%7.5f,um/min\n",				  histogram_velocity_median);	
	fprintf(VelHistFile, "Velocity mean,,%7.5f,um/min\n",				  histogram_velocity_mean  );

	// Writing DIAMETER HISTOGRAM to file 
	fprintf(DiaHistFile,"Equivalent circular diameter [um],Candidate cells,Identified cells\n");
	for( i=0; i<Nof_Bins; i++ )
		fprintf(DiaHistFile,"%7.5f,%d,%d\n", (i+0.5)*Max_Diameter_um/Nof_Bins, Diameter_Histogram[0][i], Diameter_Histogram[1][i]);

	// Calculating the averages over all valid tracks at each time step to populate the first columns of the velocity MATRIX 
	for( lev=0; lev<level; lev++ ){
		velocity_sum=velocity_sum2=0;
		m=0; 
		for( i=1; i<=nof_trax_used; i++ ){						// Loop over tracks
			if( (velocity=Velocity_Matrix[lev][5+i]) >= 0 ){
				m++;											// Counting number of valid velocities for this time step 
				velocity_sum  += velocity;						// For mean and SD calculations
				velocity_sum2 += velocity*velocity;				// For mean and SD calculations
			} 
		}
		if( m > 1 ){
			Velocity_Matrix[lev][0] =  (float)m;																	// Storing number of data points for convenience and similarity with bottom of Velocity_Matrix[][]
			Velocity_Matrix[lev][1] =  (float)sqrt( Cell_Shift_Vector_Sum_X[lev]*Cell_Shift_Vector_Sum_X[lev] + Cell_Shift_Vector_Sum_Y[lev]*Cell_Shift_Vector_Sum_Y[lev] )/(TimeStep*Image_Nr_Increment*m);			// Length of the accumulated velocity vector sum. Should be small when no image shift and large if shift correction is needed. 
			Velocity_Matrix[lev][2] =  (float)(velocity_sum/m);														// Mean velocity of all valid tracks at this time step (level)
			Velocity_Matrix[lev][3] =  (float)sqrt( (velocity_sum2-velocity_sum*velocity_sum/m)/(m-1) );			// SD  of the mean velocity at this time step (level)
			Velocity_Matrix[lev][4] =  (float)(Velocity_Matrix[lev][3]/sqrt((double)m));							// SEM of the mean velocity at this time step (level)
		} else Velocity_Matrix[lev][0]=Velocity_Matrix[lev][1]=Velocity_Matrix[lev][2]=Velocity_Matrix[lev][3]=Velocity_Matrix[lev][4]=(float)(-1);			// This should never happen, but marking as negative just in case 
	}

	// Calculating the LEVELS OVERALL AVERAGE velocity with uncertainty by averaging the time step averages over all LEVELS in the velocity matrix 
	velocity_sum=velocity_sum2=0;
	m=0; 
	for( lev=0; lev<level; lev++ ){
		if( (velocity=Velocity_Matrix[lev][2]) >= 0 ){		// Velocity_Matrix[lev][2] = average velocity at this level for all valid tracks
			m++;											// Counting number of valid values 
			velocity_sum  += velocity;						// For mean and SD calculations
			velocity_sum2 += velocity*velocity;				// For mean and SD calculations
		} 
	}
	if( m > 1 ){											// To save variables, storing overall averages over TIME STEPS (LEVELS) below data in Velocity_Matrix[][2]
		Velocity_Matrix[level+2][2] =  (float)m;																	// Mumber of valid levels
		Velocity_Matrix[level+3][2] =  (float)(velocity_sum/m);														// Mean velocity  
		Velocity_Matrix[level+4][2] =  (float)sqrt( (velocity_sum2-velocity_sum*velocity_sum/m)/(m-1) );			// SD  of velocity  
		Velocity_Matrix[level+5][2] =  (float)(Velocity_Matrix[level+4][2]/sqrt((double)m));						// SEM of velocity  
	} else Velocity_Matrix[level+2][2]=Velocity_Matrix[level+3][2]=Velocity_Matrix[level+4][2]=Velocity_Matrix[level+5][2]=(float)(-1);			// This should never happen, but marking as negative just in case 

	// Calculating the TRACKS OVERALL AVERAGE velocity with uncertainty by averaging the track averages over all TRACKS in the velocity matrix 
	velocity_sum=velocity_sum2=0;
	m=0; 
	for( i=1; i<=nof_trax_used; i++ ){
		if( (velocity=Velocity_Matrix[level+3][5+i]) >= 0 ){	// Velocity_Matrix[level+3][5+i] = average velocity for this track over all valid time steps
			m++;												// Counting number of valid values 
			velocity_sum  += velocity;							// For mean and SD calculations
			velocity_sum2 += velocity*velocity;					// For mean and SD calculations
		} 
	}
	if( m > 1 ){												// To save variables, storing overall averages over TRACKS below data in Velocity_Matrix[][1]
		Velocity_Matrix[level+2][1] =  (float)m;																// Mumber of valid levels
		Velocity_Matrix[level+3][1] =  (float)(velocity_sum/m);													// Mean velocity  
		Velocity_Matrix[level+4][1] =  (float)sqrt( (velocity_sum2-velocity_sum*velocity_sum/m)/(m-1) );		// SD  of velocity  
		Velocity_Matrix[level+5][1] =  (float)(Velocity_Matrix[level+4][1]/sqrt((double)m));						// SEM of velocity  
	} else Velocity_Matrix[level+2][1]=Velocity_Matrix[level+3][1]=Velocity_Matrix[level+4][1]=Velocity_Matrix[level+5][1]=(float)(-1);			// This should never happen, but marking as negative just in case 

	// Calculating the TRACKS OVERALL AVERAGE EUCLIDEAN velocity with uncertainty by averaging the track averages over all TRACKS in the velocity matrix 
/*	velocity_sum=velocity_sum2=0;
	m=0; 
	for( i=1; i<=nof_trax_used; i++ ){
		if( (velocity=Velocity_Matrix[level+6][5+i]) >= 0 ){	// Velocity_Matrix[level+6][5+i] = Euclidean velocity for this track over all valid time steps
			m++;												// Counting number of valid values 
			velocity_sum  += velocity;							// For mean and SD calculations
			velocity_sum2 += velocity*velocity;					// For mean and SD calculations
		} 
	} //*/
	if( m > 1 ){												// To save variables, storing overall averages over TRACKS below data in Velocity_Matrix[][1]
//		Velocity_Matrix[level+6][1] =  velocity_sum/m;													// Mean Euclidean velocity  
//		Velocity_Matrix[level+7][1] =  sqrt( (velocity_sum2-velocity_sum*velocity_sum/m)/(m-1) );		// SD  of Euclidean velocity  
//		Velocity_Matrix[level+8][1] =  Velocity_Matrix[level+7][1]/sqrt((double)m);						// SEM of Euclidean velocity  

		// Calculating the Mean, SD and SEM for Accumulated distance, Euclidean distance, Directness and the two FMI (Forward Migration Index) over all tracks, to be reported in _Overall_summary.csv 
		accumulated_distance_mean = accumulated_distance_sum/m;
		accumulated_distance_sd   = sqrt( (accumulated_distance_sum2-accumulated_distance_sum*accumulated_distance_sum/m)/(m-1) );
		accumulated_distance_sem  = accumulated_distance_sd/sqrt((double)m);

		euclidean_distance_mean = euclidean_distance_sum/m;
		euclidean_distance_sd   = sqrt( (euclidean_distance_sum2-euclidean_distance_sum*euclidean_distance_sum/m)/(m-1) );
		euclidean_distance_sem  = euclidean_distance_sd/sqrt((double)m);

		euclidean_velocity_mean = euclidean_velocity_sum/m;
		euclidean_velocity_sd   = sqrt( (euclidean_velocity_sum2-euclidean_velocity_sum*euclidean_velocity_sum/m)/(m-1) );
		euclidean_velocity_sem  = euclidean_velocity_sd/sqrt((double)m);

		directness_mean = directness_sum/m;
		directness_sd   = sqrt( (direcntess_sum2-directness_sum*directness_sum/m)/(m-1) );
		directness_sem  = directness_sd/sqrt((double)m);

		fmi_x_mean = fmi_x_sum/m;
		fmi_x_sd   = sqrt( (fmi_x_sum2-fmi_x_sum*fmi_x_sum/m)/(m-1) );
		fmi_x_sem  = fmi_x_sd/sqrt((double)m);

		fmi_y_mean = fmi_y_sum/m;
		fmi_y_sd   = sqrt( (fmi_y_sum2-fmi_y_sum*fmi_y_sum/m)/(m-1) );
		fmi_y_sem  = fmi_y_sd/sqrt((double)m);

	} else accumulated_distance_mean=accumulated_distance_sd=accumulated_distance_sem=euclidean_distance_mean=euclidean_distance_sd=euclidean_distance_sem=euclidean_velocity_mean=euclidean_velocity_sd=euclidean_velocity_sem=directness_mean=directness_sd=directness_sem=fmi_x_mean=fmi_x_sd=fmi_x_sem=fmi_y_mean=fmi_y_sd=fmi_y_sem=-1;		// This should never happen, but marking as negative just in case 

	// Calculating the TOTAL OVERALL AVERAGE velocity with uncertainty by averaging separately ALL CELL MOVEMENTS in valid tracks from the velocity matrix 
	velocity_sum=velocity_sum2=0;
	m=0; 	
	for( lev=0; lev<level; lev++ ){
		for( i=1; i<=nof_trax_used; i++ ){
			if( (velocity=Velocity_Matrix[lev][5+i]) >= 0 ){	// Velocity_Matrix[lev][5+i] = velocity of this particular cell movement in a valid track
				m++;											// Counting number of valid values 
				velocity_sum  += velocity;						// For mean and SD calculations
				velocity_sum2 += velocity*velocity;				// For mean and SD calculations
			} 
		}
	}
	if( m > 1 ){												// To save variables, storing overall averages over MOVEMENTS below data in Velocity_Matrix[][0]
		Velocity_Matrix[level+2][0] =  (float)m;																// Mumber of valid levels
		Velocity_Matrix[level+3][0] =  (float)(velocity_sum/m);													// Mean velocity  
		Velocity_Matrix[level+4][0] =  (float)sqrt( (velocity_sum2-velocity_sum*velocity_sum/m)/(m-1) );		// SD  of velocity  
		Velocity_Matrix[level+5][0] =  (float)(Velocity_Matrix[level+4][0]/sqrt((double)m));					// SEM of velocity  
	} else Velocity_Matrix[level+2][0]=Velocity_Matrix[level+3][0]=Velocity_Matrix[level+4][0]=Velocity_Matrix[level+5][0]=(float)(-1);			// This should never happen, but marking as negative just in case 

	// Writing velocity MATRIX to file 
	fprintf(VelMatrFile, ",,,,Averaging velocities over,,,First image,");																															for( i=1; i<=nof_trax_used; i++ ) fprintf(VelMatrFile,"%1.0f,", Velocity_Matrix[level  ][5+i]); fprintf(VelMatrFile,"\n"); 
	fprintf(VelMatrFile, ",,,Movements,Cell Tracks,Time Steps,,Last image,");																														for( i=1; i<=nof_trax_used; i++ ) fprintf(VelMatrFile,"%1.0f,", Velocity_Matrix[level+1][5+i]); fprintf(VelMatrFile,"\n"); 
	fprintf(VelMatrFile, ",,Number of values,%7.5f,%7.5f,%7.5f,,Number of cell movements,",			  Velocity_Matrix[level+2][0], Velocity_Matrix[level+2][1], Velocity_Matrix[level+2][2]);		for( i=1; i<=nof_trax_used; i++ ) fprintf(VelMatrFile,"%1.0f,", Velocity_Matrix[level+2][5+i]); fprintf(VelMatrFile,"\n"); 
	fprintf(VelMatrFile, ",Overall velocity mean,[um/min],%7.5f,%7.5f,%7.5f,Velocity Mean,[um/min],", Velocity_Matrix[level+3][0], Velocity_Matrix[level+3][1], Velocity_Matrix[level+3][2]);		for( i=1; i<=nof_trax_used; i++ ) fprintf(VelMatrFile,"%7.5f,", Velocity_Matrix[level+3][5+i]); fprintf(VelMatrFile,"\n"); 
	fprintf(VelMatrFile, ",Overall velocity   SD,[um/min],%7.5f,%7.5f,%7.5f,Velocity   SD,[um/min],", Velocity_Matrix[level+4][0], Velocity_Matrix[level+4][1], Velocity_Matrix[level+4][2]);		for( i=1; i<=nof_trax_used; i++ ) fprintf(VelMatrFile,"%7.5f,", Velocity_Matrix[level+4][5+i]); fprintf(VelMatrFile,"\n"); 
	fprintf(VelMatrFile, ",Overall velocity  SEM,[um/min],%7.5f,%7.5f,%7.5f,Velocity  SEM,[um/min],", Velocity_Matrix[level+5][0], Velocity_Matrix[level+5][1], Velocity_Matrix[level+5][2]);		for( i=1; i<=nof_trax_used; i++ ) fprintf(VelMatrFile,"%7.5f,", Velocity_Matrix[level+5][5+i]); fprintf(VelMatrFile,"\n"); 

	fprintf(VelMatrFile, "\n");
	fprintf(VelMatrFile, "Image name,Time,Time,Number of,Velocity sum / cell,Velocity mean,Velocity SD,Velocity SEM,");		// Headers for summary data averaged over each time step
	for( i=1; i<=nof_trax_used; i++ ) fprintf(VelMatrFile,"Track #,"); 
	if( Wound_Healing_Mode ) fprintf(VelMatrFile, "Mean overall velocity sum X,Mean overall velocity sum Y,Mean upper half velocity sum X,Mean upper half velocity sum Y,Mean lower half velocity sum X,Mean lower half velocity sum Y,Mean velocity into wound,Number of tracks,Number of tracks\n");										// Units for summary data averaged over each time step
	else fprintf(VelMatrFile, "\n");	 
	fprintf(VelMatrFile, ",[min],[h],tracks,[um/min],[um/min],[um/min],[um/min],");										// Units for summary data averaged over each time step
	for( i=1; i<=nof_trax_used; i++ ) fprintf(VelMatrFile,"%d,", i); 
	if( Wound_Healing_Mode ) fprintf(VelMatrFile, "[um/min],[um/min],[um/min],[um/min],[um/min],[um/min],[um/min],in upper half,in lower half\n");
	else fprintf(VelMatrFile, "\n");	 
	for( lev=0; lev<level-1; lev++ ){									// One less velocity than images since calculated from one image to the one after 
		fprintf(VelMatrFile, "%s%04d,%3.1f,%4.2f,", BMPfilename, lev+First_Image_Number, (lev+First_Image_Number)*TimeStep*Image_Nr_Increment, (lev+First_Image_Number)*TimeStep*Image_Nr_Increment/60.0);
		for( i=0; i<5;				i++ )																   fprintf(VelMatrFile,"%7.5f,", Velocity_Matrix[lev][  i] );		// First writing summary statistics which is stored at right edge of matrix. Not checking for negative numbers here since Velocity sum X and Y may be negative.
		for( i=1; i<=nof_trax_used; i++ ) if( Velocity_Matrix[lev][5+i]<0 ) fprintf(VelMatrFile,","); else fprintf(VelMatrFile,"%7.5f,", Velocity_Matrix[lev][5+i] );		// Then writing matrix data for more convenient processing in Excel with summary data first in rows and columns 
		if( Wound_Healing_Mode ) fprintf(VelMatrFile,"%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%d,%d\n", Cell_Shift_Vector_Sum_X[lev]/(TimeStep*Image_Nr_Increment*Velocity_Matrix[lev][0]), Cell_Shift_Vector_Sum_Y[lev]/(TimeStep*Image_Nr_Increment*Velocity_Matrix[lev][0]), Cell_Shift_Vector_Sum_X_WH_Upper[lev]/(TimeStep*Image_Nr_Increment*WH_Nof_Cells_Upper[lev]), Cell_Shift_Vector_Sum_Y_WH_Upper[lev]/(TimeStep*Image_Nr_Increment*WH_Nof_Cells_Upper[lev]), Cell_Shift_Vector_Sum_X_WH_Lower[lev]/(TimeStep*Image_Nr_Increment*WH_Nof_Cells_Lower[lev]), Cell_Shift_Vector_Sum_Y_WH_Lower[lev]/(TimeStep*Image_Nr_Increment*WH_Nof_Cells_Lower[lev]), (-Cell_Shift_Vector_Sum_Y_WH_Upper[lev]+Cell_Shift_Vector_Sum_Y_WH_Lower[lev])/(TimeStep*Image_Nr_Increment*(WH_Nof_Cells_Upper[lev]+WH_Nof_Cells_Lower[lev])), WH_Nof_Cells_Upper[lev], WH_Nof_Cells_Lower[lev]); 
		else					 fprintf(VelMatrFile,"\n"); 
	}

	// Finding well coordinate letter and number from BMPfilename
	j=l=0; 
	while(l<strlen(BMPfilename) && (ch=BMPfilename[l])!='_') {
		if( ch>='0' && ch<='9' ) video_number_string[j++]=ch;								// Copying all digits in image name prior to first _ over to video_number_string for later conversion to integer
		l++;
	}
	sscanf(video_number_string, "%d", &video_number);										// Converting string of digits to integer
	well_letter = BMPfilename[l+1];															// Reading well coordinate letter as the character after the first underscore
	well_number = BMPfilename[l+2]-'0';														// Taking the well coordinate number as the digit following immediately after 
	if( (ch=BMPfilename[l+3])>='0' && ch<='9' ) well_number = 10*well_number + ch-'0';		// Checking if number has two digits and then reading it correctly 

	// Summary data for entire video written as one line in the _Overall_summary.csv file
	fprintf(SummaryFile,"%s%04d,%d,%s%04d,%d,%1.0f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%1.0f,%7.5f,%7.5f,%7.5f,%1.0f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%7.5f,%d,%c,%d\n", BMPfilename, First_Image_Number, Nof__Cells[0], BMPfilename, First_Image_Number+level-1, Nof__Cells[level-1],  Velocity_Matrix[level+2][1], Velocity_Matrix[level+3][1], Velocity_Matrix[level+4][1], Velocity_Matrix[level+5][1], directness_mean, directness_sd, directness_sem, fmi_x_mean,  fmi_x_sd, fmi_x_sem, fmi_y_mean, fmi_y_sd, fmi_y_sem, accumulated_distance_mean, accumulated_distance_sd, accumulated_distance_sem, euclidean_distance_mean, euclidean_distance_sd, euclidean_distance_sem, euclidean_velocity_mean, euclidean_velocity_sd, euclidean_velocity_sem, Velocity_Matrix[level+2][2], Velocity_Matrix[level+3][2], Velocity_Matrix[level+4][2], Velocity_Matrix[level+5][2], Velocity_Matrix[level+2][0], Velocity_Matrix[level+3][0], Velocity_Matrix[level+4][0], Velocity_Matrix[level+5][0], histogram_velocity_mean, histogram_velocity_mode_parabola, histogram_velocity_median, histogram_velocity_75_percentile, video_number, well_letter, well_number);	
	fprintf(CellCountFile,"%s,%d,%c,%d,# cells in image,", BMPfilename, video_number, well_letter, well_number); 
	for( n=0; n<First_Image_Number;          n++ ) fprintf(CellCountFile, ",");																		// Writing empty columns if this video starts
	for( n=First_Image_Number; n<=First_Image_Number+level-1; n++ ) fprintf(CellCountFile, "%d,", Nof_Identified_Cells[n-First_Image_Number]);		// In the rare cases where First_Image_Number+level-1 > Tuning__Last_Image, there will not be time values above the last columns of counts for this row. But better to give all cell counts than to trunkate.
	fprintf(CellCountFile, "\n"); 

	fclose(PosFile);
	fclose(MatrixFile);
	fclose(TraxFile);
	fclose(DiaHistFile);
	fclose(VelHistFile);
	fclose(VelMatrFile);
	fclose(ShiftFile);
}


// Reading back shifted images at each level, drawing tracks and color dots in all cells with long enough tracks 

void generate_valid_tracks_images(int level)
{
	char filename[STRING_LENGTH];
	unsigned char color[3]; 
	int i, k, l, m, lev; 
	int int_radius, x, y, x_1, y_1, x_2, y_2, m_1, m_2, d, delta_x, delta_y, delta_max;	
	float dot_radius;
	double ssum, smean, greyvalue; 

	for( lev=0; lev<level; lev++ ){

	printf("Generating final video image %3d of %3d...     \r", lev, level-1); 

		sprintf(filename, "%s/%s_01_%04d_Shifted.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*lev);
		read_bitmap(Oimage, filename);
		if( !Write_1_SH_img ) remove(filename);				// Deleting shifted images after use if command file does not call for them to be stored 

		// Simple smoothing filter to make images with high contrast regions look more similar if some are shifted much and thus pixels are mixed more in the fractional shifting above. 
		ssum = 0; 
		for( y=1; y<Bimage->ysize-1; y++ )
			for( x=1; x<Bimage->xsize-1; x++ ) {
				greyvalue = 0.4*Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][0] + 0.15*Oimage->image[Left___Crop_Margin+x+1][Bottom_Crop_Margin+y][0] + 0.15*Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y+1][0] + 0.15*Oimage->image[Left___Crop_Margin+x-1][Bottom_Crop_Margin+y][0] + 0.15*Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y-1][0]; 
				ssum  += greyvalue;
				Smoothed_Image[Left___Crop_Margin+x][Bottom_Crop_Margin+y] = (float)greyvalue;					// Using Smoothed_Image[][] as temporary storage to keep grey levels as floats with no roundoff to integers 
			}  //*/
		smean = ssum/(Bimage->xsize*Bimage->ysize);									// Should never give division by zero error so not checking  

		// Calculating image with option to increase contrast through the factor VT_Image_Contrast while keeping the mean grey level as it was 
		for( y=1; y<Bimage->ysize-1; y++ )
			for( x=1; x<Bimage->xsize-1; x++ ) {
				greyvalue = smean + VT_Image_Contrast*( Smoothed_Image[Left___Crop_Margin+x][Bottom_Crop_Margin+y] - smean ); 
				if(      greyvalue <   0 ) greyvalue =   0;
				else if( greyvalue > 255 ) greyvalue = 255;
				for( k=0; k<3; k++) Oimage->image[Left___Crop_Margin+x][Bottom_Crop_Margin+y][k] = (unsigned char)(greyvalue + 0.5);  
			}  		

		// Generating image overlays with trace lines and color dots in all cells with valid tracks
		convert2color(Oimage);																			// To allow showing different cells with different colors 
		for( i=1; i<=Nof_Matched_Cells; i++ ){

			if( Matched_Cells[level][i] >= Shortest_Cell_Track ){										// Only processing tracks that have been matched for at least Shortest_Cell_Track images 

//				printf("Matched cell number %d has %d images in the track\n", i, Matched_Cells[level][i]); 
				m = Cells[lev][Matched_Cells[lev][i]].cellnumber;										// The cell number used to get the same color as the cell in the image 
				cellnumber2color(m,color); 

				// Drawing the trace of all valid tracks with the line of the same color as the cell, to see more easily how much the cell has moved 
				if( Draw_Track && lev > 0 ){
					for( l=1; l<=lev; l++ ){
						if( m>0 && (m_2=Matched_Cells[l][i])>0 && (m_1=Matched_Cells[l-1][i])>0 ){		// We have a matched cell in this lev and the previous level
							x_2 = Cells[l  ][m_2].cx;													// Shorthand 
							y_2 = Cells[l  ][m_2].cy;													// Unit: [pixels]
							x_1 = Cells[l-1][m_1].cx;			
							y_1 = Cells[l-1][m_1].cy;
							delta_x = x_2 - x_1; 
							delta_y = y_2 - y_1; 
							delta_max = maxi(abs(delta_x), abs(delta_y));
							if( delta_max == 0 ) delta_max = 1; 
							for( d=0; d<=delta_max; d++ ) { 
								for( k=0; k<3; k++) Oimage->image[(int)(Left___Crop_Margin+x_1+1.0*d*delta_x/delta_max+0.5)][(int)(Bottom_Crop_Margin+y_1+1.0*d*delta_y/delta_max+0.5)][k] = color[k];
							}
						}
					}
				}

				// Drawing a filled circle with the correct cell color at the center of each cell in the original image 
				dot_radius = (float)(0.5*Dot_Diameter_um/Pixel_Size);		// Radius in um  
				int_radius = (int)(dot_radius+0.5);							// Rounding off the diameter to find integer limits to the two loops below  
				if( int_radius < 1 ) int_radius = 1;						// Using at least 3 (=1+2*int_radius) as smallest spot diameter to draw in original image at cell center of mass
				if( Draw_CM && Matched_Cells[lev][i]>0 ){					// Only drawing dot at those levels that were sucessfully  matched 
					for( y=-int_radius; y<=int_radius; y++ )
						for( x=-int_radius; x<=int_radius; x++ )
							if( sqrt(1.0*x*x+y*y)<dot_radius && Left___Crop_Margin+Cells[lev][i].cx+x>0 && Bottom_Crop_Margin+Cells[lev][i].cy+y>0 && Left___Crop_Margin+Cells[lev][i].cx+x<Oimage->xsize && Bottom_Crop_Margin+Cells[lev][i].cy+y<Oimage->ysize ) {
								for( k=0; k<3; k++) Oimage->image[Left___Crop_Margin+Cells[lev][Matched_Cells[lev][i]].cx+x][Bottom_Crop_Margin+Cells[lev][Matched_Cells[lev][i]].cy+y][k] = color[k]; 	
							}
//					for( k=0; k<3; k++) Oimage->image[(int)(Left___Crop_Margin+Cells[lev][Matched_Cells[lev][i]].fx+0.5)][(int)(Bottom_Crop_Margin+Cells[lev][Matched_Cells[lev][i]].fy+0.5)][k] = (unsigned char)0; 	// Black dot at original center before time filtering just for comparison
				}
			}
		}	// End of for (i...) over all matched cells 

		// Writing image with only VALID TRACKS included for making video or for checking, including scalebar if desired
		if( Add_Scalebar ) draw_scalebar(); 
		sprintf(filename, "%s/%s_07_%04d_Valid_Track.bmp", ResultsFolderName, BMPfilename, First_Image_Number+Image_Nr_Increment*lev);
		out_bitmap(Oimage, filename);

	}  // End of for( lev...) over all levels 
	printf("Generated %d images where only valid tracks are shown.                                          \n", level); 
}


// Setting all arrays to zero, thereby removing all memory of the previous set of images 

void reset_arrays(void){

	int i, j; 

	for( i=0; i<2; i++ )
		for( j=0; j<MAXBINS; j++ )
			Diameter_Histogram[i][j]=0; 

	for( i=0; i<MAXLEVELS; i++ )
		for( j=0; j<MAXBINS; j++ )
			Velocity_Histogram[i][j]=0; 

	for( i=0; i<MAXLEVELS; i++ )
		for( j=0; j<MATCHED_CELLS; j++ ){
			Matched_Cells[i][j]=0;
			Velocity_Matrix[i][j]=0; 
		}

	for( i=0; i<MAXLEVELS; i++ ){
		Nof__Cells[i]=0; 
		Image_Shift_um_X[i]=0;
		Image_Shift_um_Y[i]=0;
		Cell_Shift_Vector_Sum_X[i]=0;
		Cell_Shift_Vector_Sum_Y[i]=0;
		Cell_Shift_Vector_Sum_X_WH_Upper[i]=0;
		Cell_Shift_Vector_Sum_Y_WH_Upper[i]=0;
		Cell_Shift_Vector_Sum_X_WH_Lower[i]=0;
		Cell_Shift_Vector_Sum_Y_WH_Lower[i]=0;
		WH_Nof_Cells_Upper[i]=0;
		WH_Nof_Cells_Lower[i]=0;
	}

	for( i=0; i<MAXCELLS+1; i++ )
		for( j=0; j<MAXCELLS; j++ )
			Matching_Matrix[i][j]=0; 

	for( i=0; i<MAXLEVELS; i++ )
		for( j=0; j<MAXCELLS; j++ ){
			Cells[i][j].area_fraction	=0; 
			Cells[i][j].area_px			=0; 
			Cells[i][j].aspect_ratio	=0; 
			Cells[i][j].cellnumber		=0; 
			Cells[i][j].cm_x_first		=0; 
			Cells[i][j].cm_y_first		=0; 
			Cells[i][j].cx				=0; 
			Cells[i][j].cy				=0; 
			Cells[i][j].dx				=0; 
			Cells[i][j].dy				=0; 
			Cells[i][j].max_x			=0; 
			Cells[i][j].max_y			=0; 
			Cells[i][j].min_x			=0; 
			Cells[i][j].min_y			=0; 
		}
	
	First_Image_Number = Last__Image_Number = Nof_Matched_Cells = 0;		// Also resetting these series specific variables 
	Previous_Image_Shift_X_px = Previous_Image_Shift_Y_px = 0;

	fprintf(CMDWindowFile, "\n");	// Just putting one blank line after each video in the .csv file to get visible grouping and easier plotting in Excel 
}


// Drawing a 100 µm scalebar and text in output image if the two flags are true

void draw_scalebar(void)
{
	unsigned char paste_color[3];
	int scalebar_length_px, scalebar_height_px, h_um_image_length_px, corner_margins_px, scalebar_and_image_midpoint_x_px;
	int k, x, y, startx, starty;	
	long g, gg[7]={936350, 608848, 608848, 948816, 608848, 608848, 936336}; 

	h_um_image_length_px = 100; 
	corner_margins_px = 10;
	if( Add_Scalebar ){
		translate_color(Scalebar_Color, paste_color);
		scalebar_length_px = (int)(100/Pixel_Size     + 0.5);
		scalebar_height_px = (int)(scalebar_length_px/10 + 0.5);
		scalebar_and_image_midpoint_x_px = (int)(Oimage->xsize - corner_margins_px - 0.5*maxi(h_um_image_length_px, scalebar_length_px) + 0.5);	// Using the longer of the scale bar or "100 µm" image to find the common midpoint 
		startx = (int)(scalebar_and_image_midpoint_x_px - 0.5*scalebar_length_px + 0.5);			// Centering scalebar and "100 µm" image horizontally
		starty = corner_margins_px;																	// Lower right corner of scale bar lies corner_margins_px from image corner in both x and y
		for( y=0; y<scalebar_height_px; y++ )
			for( x=0; x<scalebar_length_px; x++ )
				for( k=0; k<3; k++) Oimage->image[startx+x][starty+y][k] = paste_color[k];
	}
	for( y=0; y<7; y++ ){	
		g = gg[y]; x = 19;
		while( (g=g>>1)>0 ) { if( g&1 && Oimage->image[x][y+1][0]<255 ) Oimage->image[x][y+1][0]+=1; x--; }
	}
	if( Add_100um_image ){
		startx = (int)(scalebar_and_image_midpoint_x_px - 0.5*Simage->xsize + 0.5);					// Centering scalebar and "100 µm" image horizontally
		starty = corner_margins_px + scalebar_height_px;											// Starting "100 µm" image just above scalebar line
		paste_grayscale_with_white_transparent(Simage, Oimage, paste_color, startx, starty);
	}
}


// Writing error message to screen and to file, quitting if not a warning 
void write_error_message(FILE *fp, char message_type, char *message)
{
	if( message_type == WARNING ) { printf("Warning in CellTraxx: %s\n", message);		fprintf(fp, "Warning in celltraxx.exe: %s\n", message); }
	else				          { printf("\nError in CellTraxx: %s\n", message);		fprintf(fp, "Error in celltraxx.exe: %s",     message); exit(0); }
}


// Smooths gray pixels in the bitmap with a gaussian filter of size filterradius (= fraction of ysize). 

int gaussian_populate_filter(int filterradius)
{
	char message[STRING_LENGTH];
	int fx, fy, gauss_halfsize;
	float k, filtercutoff, filtersum;

	if( filterradius<2  )			{ sprintf(message, "The variable filterradius is less than 2.\n");																						write_error_message(ErrorMessageFile, ERROR, message); }
	if( filterradius>0.25*MAXY )	{ sprintf(message, "The variable filterradius %d is too large. Max value is a quarter of the image height (%3.0f pixels).\n", filterradius, 0.25*MAXY);	write_error_message(ErrorMessageFile, ERROR, message); }
	if( filterradius>=SHIFTXY)		{ sprintf(message, "The variable filterradius %d is too large. It must be smaller than the image margin (%d pixels).\n",		filterradius, SHIFTXY);	write_error_message(ErrorMessageFile, ERROR, message); }
	if( filterradius>=MAXX-Bimage->xsize-SHIFTXY  || filterradius>=MAXY-Bimage->ysize-SHIFTXY )		{ sprintf(message, "The variable filterradius %d is larger than the top image margin (%d pixels) or right image margin (%d pixels).\n", filterradius, MAXX-Bimage->xsize-SHIFTXY, MAXY-Bimage->ysize-SHIFTXY );  write_error_message(ErrorMessageFile, ERROR, message); }

	filtercutoff = (float)0.01;														// Ending Gaussian when it reaches 0.01 
	gauss_halfsize = filterradius-1;												// in pixels 
	if(gauss_halfsize<2) gauss_halfsize=2;		
	k=(float)(sqrt(-log(filtercutoff))/filterradius);								// scaling factor 

	filtersum=0;
	for(fx=-gauss_halfsize; fx<=gauss_halfsize; fx++)								// Calculating Gaussian filter 
		for(fy=-gauss_halfsize; fy<=gauss_halfsize; fy++)
			filtersum += ( Gaussian_Filter[fx+gauss_halfsize][fy+gauss_halfsize]=(float)exp(-fx*fx*k*k-fy*fy*k*k) );

	for(fx=-gauss_halfsize; fx<=gauss_halfsize; fx++)								// Normalizing filter 
		for(fy=-gauss_halfsize; fy<=gauss_halfsize; fy++)
			Gaussian_Filter[fx+gauss_halfsize][fy+gauss_halfsize] /= filtersum;

	return( gauss_halfsize );
}


// Smooths the Expanded_Image[][] by the Gaussian_Filter[][] over to the Smoothed_Image[][] 

void gaussian_smooth(int gauss_halfsize)
{
	int x, y, xs,  ys;
	int fx, fy;
	float meanvalue;

	xs=SHIFTXY+Bimage->xsize; 		// Shorthand to save calculation time
	ys=SHIFTXY+Bimage->ysize; 

	for(x=SHIFTXY; x<xs; x++){
		for(y=SHIFTXY; y<ys; y++){
			meanvalue=0;
			for(fx=-gauss_halfsize; fx<=gauss_halfsize; fx++)
				for(fy=-gauss_halfsize; fy<=gauss_halfsize; fy++)
					meanvalue += Expanded_Image[x+fx][y+fy]*Gaussian_Filter[gauss_halfsize+fx][gauss_halfsize+fy];	
			Smoothed_Image[x][y]=meanvalue;
		}
	}
}



// Flood fills Segmented_Image[][] with the given number as long as four-connected pixels of the same number are available. Returns the number of pixels filled, and also the four bounding box values by reference.

int fill_segmented_box(int fillnumber, int withnumber, int startx, int starty, int *x_min, int *x_max, int *y_min, int *y_max)  
{
	int x, y;
	int plistread, plistwrite, count;

	plistread=plistwrite=count=0;

	*x_min = *x_max = startx;					// Initially setting the bounding box values equal to the starting point
	*y_min = *y_max = starty;

	FILL_SEG(startx, starty)					// Macro for FILL_SEG is #defined at the top 

	while( plistread!=plistwrite ){
		x=PointList_x[plistread];
		y=PointList_y[plistread];
		plistread=(plistread+1)%(2*MAXY);

		if( x < *x_min ) *x_min = x;			// Checking if current point expands the bounding box 
		if( x > *x_max ) *x_max = x; 
		if( y < *y_min ) *y_min = y; 
		if( y > *y_max ) *y_max = y; 

		FILL_SEG(x+1, y  )
		FILL_SEG(x  , y-1)
		FILL_SEG(x-1, y  )
		FILL_SEG(x  , y+1)
	}
	return(count);
}


// Flood fills Segmented_Image[][] with the given number as long as four-connected pixels of the same number are available. Returns the number of pixels filled. 

int fill_segmented(int fillnumber, int withnumber, int startx, int starty)  
{
	int x, y;
	int plistread, plistwrite, count;

	plistread=plistwrite=count=0;

	FILL_SEG(startx, starty)					// Macro for FILL_SEG is #defined at the top 

	while( plistread!=plistwrite ){
		x=PointList_x[plistread];
		y=PointList_y[plistread];
		plistread=(plistread+1)%(2*MAXY);

		FILL_SEG(x+1, y  )
		FILL_SEG(x  , y-1)
		FILL_SEG(x-1, y  )
		FILL_SEG(x  , y+1)
	}
	return(count);
}


// Fills Tuning_Image[][] with the given number as long as four-connected pixels with fillnumber are available. Returns the number of pixels filled.

int fill_tuning(int fillnumber, int withnumber, int startx, int starty)  
{
	int x, y;
	int plistread, plistwrite, count;

	plistread=plistwrite=count=0;

	FILL_TUN(startx, starty)					// Macro for FILL_TUN is #defined at the top 

	while( plistread!=plistwrite ){
		x=PointList_x[plistread];
		y=PointList_y[plistread];
		plistread=(plistread+1)%(2*MAXY);

		FILL_TUN(x+1, y  )
		FILL_TUN(x  , y-1)
		FILL_TUN(x-1, y  )
		FILL_TUN(x  , y+1)
	}
	return(count);
}


// Returns a color value for each channel which differs significantly from the color returned by the neighboursing cellnumbers 

 void cellnumber2color(short cellnumber,  unsigned char color[3])
{
	unsigned char blue, green, red;
	double hue, lum, b, g, r;

	if( cellnumber < 0 ) cellnumber = -cellnumber; 

	hue = (float)( (cellnumber*137)%(6*HUESPERCOLOR) );							// Picks a different hue for each neighboursing cellnumber between 0 and 120 (=6*HUESPERCOLOR) from which to calculate the hue along the spectrum blue - cyan - green - yellow - red - magenta - blue
	if( Color_Print_Out_Mode ) lum = (float)( 220 + 35*(cellnumber&3)/3.0 );	// Will take values 220, 232, 243 or 255 to give various degree of darkness for each hue
	else                       lum = (float)( ((cellnumber&3)+5)*255.0/8 );		// Will take values 159, 191, 223 or 255 to give various degree of darkness for each hue

	b = (absd(3 - hue/HUESPERCOLOR) - 1);	if( b < 0 ) b = 0; else if( b > 1 ) b = 1;	blue  = (unsigned char)(b*lum);
	g = (2 - absd(2 - hue/HUESPERCOLOR)); 	if( g < 0 ) g = 0; else if( g > 1 ) g = 1;	green = (unsigned char)(g*lum);
	r = (2 - absd(4 - hue/HUESPERCOLOR)); 	if( r < 0 ) r = 0; else if( r > 1 ) r = 1;	red   = (unsigned char)(r*lum);

	// Some special cases treated as exceptions 
	     if( cellnumber == 1 ) blue = green = red = (unsigned char)(  1);		// Just to get black regions in segmented and eroded image  
	else if( cellnumber == 0 ) blue = green = red = (unsigned char)(255);		// Making background white 
	else if( cellnumber <  0 ) blue = green = red = (unsigned char)(  1);		// Making non-matched cells as black 

	color[0] = blue;
	color[1] = green;
	color[2] = red;
}


// Returns a monochrome color value with different luminocities depending on cellnumber. 

void cellnumber2monocolor(short cellnumber, unsigned char color[3], char colorcode)
{
	double lum, c, m;

	if( cellnumber < 0 ) cellnumber = -cellnumber; 

	if( Color_Print_Out_Mode ) lum = 0.8 + 0.2*(cellnumber&7)/7;		// Luminosity values from 80 % to 100 % in 8 different steps. 
	else                       lum = 0.4 + 0.4*(cellnumber&7)/7;		// Luminosity values from 40 % to  80 % in 8 different steps. 
	c = 1 - absd(2*lum-1);
	m = lum - c/2; 

	color[0] = (unsigned char)(255*(m+c*(colorcode&1)/1));		// Blue	
	color[1] = (unsigned char)(255*(m+c*(colorcode&2)/2));		// Green				colorcode: Blue = bit 1, Green = bit 2, Red = bit 3
	color[2] = (unsigned char)(255*(m+c*(colorcode&4)/4));		// Red 
}


// Reads a bitmap directly from a file assumed to exist, returns FALSE if something goes wrong	

extern char read_bitmap(BITMAP *bimp, char *infilename)
{
	unsigned char dummy, colornumber;
	char message[STRING_LENGTH]; 
	int x, y, i, xzeros;
	FILE *fp;												// Pointer to the input file 

	if( !splitbmpname(infilename, bimp) ){ sprintf(message, "Invalid path '%s', file name '%s' or extension '%s'.\n", bimp->path, bimp->filename, bimp->ext); write_error_message(ErrorMessageFile, ERROR, message); }

	sprintf(infilename, "%s%s.%s", bimp->path, bimp->filename, bimp->ext); 
	if( strcmp(bimp->ext,imgext) ) { sprintf(message, "Wrong extension (%s differs from %s).\n", bimp->ext, imgext); write_error_message(ErrorMessageFile, ERROR, message);  }

	if( (fp=fopen(infilename, "rb"))==NULL ) { sprintf(message, "The file '%s' does not exist.\n", infilename); write_error_message(ErrorMessageFile, ERROR, message); }

	bimp->b=getc(fp);			// reading signature		
	bimp->m=getc(fp);													
	if( bimp->b!='B' && bimp->m!='M' ) { sprintf(message, "Invalid bitmap signature '%c%c'.\n", bimp->b, bimp->m); write_error_message(ErrorMessageFile, ERROR, message); }
	if( fread(&bimp->filesize,       4, 1, fp)<1 ) write_error_message(ErrorMessageFile, ERROR, "Can't read file size.\n");	 
	if( fread(&bimp->reserved,       4, 1, fp)<1 ) write_error_message(ErrorMessageFile, ERROR, "Can't read reserved.\n");								 
	if( fread(&bimp->dataoffset,     4, 1, fp)<1 ) write_error_message(ErrorMessageFile, ERROR, "Can't read dataoffset.\n");								 
	if( fread(&bimp->infoheadersize, 4, 1, fp)<1 ) write_error_message(ErrorMessageFile, ERROR, "Can't read infoheadersize.\n");								 
	if( fread(&bimp->xsize,			 4, 1, fp)<1 ) write_error_message(ErrorMessageFile, ERROR, "Can't read xsize.\n");								 
	if( fread(&bimp->ysize,		     4, 1, fp)<1 ) write_error_message(ErrorMessageFile, ERROR, "Can't read ysize.\n");								 
	if( fread(&bimp->planes,		 2, 1, fp)<1 ) write_error_message(ErrorMessageFile, ERROR, "Can't read planes.\n");								 
	if( fread(&bimp->bitcount,		 2, 1, fp)<1 ) write_error_message(ErrorMessageFile, ERROR, "Can't read bitcount.\n");								 
	if( fread(&bimp->compression,	 4, 1, fp)<1 ) write_error_message(ErrorMessageFile, ERROR, "Can't read compression.\n");								 
	if( fread(&bimp->imagesize,      4, 1, fp)<1 ) write_error_message(ErrorMessageFile, ERROR, "Can't read imagesize.\n");								 
	if( fread(&bimp->xres,		     4, 1, fp)<1 ) write_error_message(ErrorMessageFile, ERROR, "Can't read xres.\n");								 
	if( fread(&bimp->yres,		     4, 1, fp)<1 ) write_error_message(ErrorMessageFile, ERROR, "Can't read yres.\n");								 
	if( fread(&bimp->colorsused,     4, 1, fp)<1 ) write_error_message(ErrorMessageFile, ERROR, "Can't read colorsused.\n");								 
	if( fread(&bimp->colorsimport,   4, 1, fp)<1 ) write_error_message(ErrorMessageFile, ERROR, "Can't read colorsimport.\n");								 

	if( bimp->xsize > MAXX ){ sprintf(message, "Image x-size %d exceedes upper limit %d.\nIncrease MAXX and recompile.\n", bimp->xsize, MAXX); write_error_message(ErrorMessageFile, ERROR, message); }								 
	if( bimp->ysize > MAXY ){ sprintf(message, "Image y-size %d exceedes upper limit %d.\nIncrease MAXY and recompile.\n", bimp->ysize, MAXY); write_error_message(ErrorMessageFile, ERROR, message); }								 

	if( bimp->bitcount==24 && bimp->compression==0 ){		// It's an uncompressed color image 
		xzeros = (bimp->xsize)%4;							// Finding number of zeros buffered at end of each row. All rows contain 4n bytes 
		for(y=0; y<bimp->ysize; y++){
			for(x=0; x<bimp->xsize; x++) if( fread(&bimp->image[x][y][0], 1, 3, fp)<3 ){ sprintf(message, "Can't read pixel [%d][%d].\n", x+1, y+1); write_error_message(ErrorMessageFile, ERROR, message); }		
			for(i=0; i<xzeros; i++)      if( fread(&dummy, 1, 1, fp               )<1 ){ sprintf(message, "Can't read zeros after line %d.\n", y+1); write_error_message(ErrorMessageFile, ERROR, message); }
		}
	} else if( bimp->bitcount==8 && bimp->compression==0 ){	// It's an uncompressed gray scale image 
		xzeros = 4-(bimp->xsize)%4;	
		if( xzeros==4 ) xzeros=0;							// Finding number of zeros buffered at end of each row. All rows contain 4n bytes 
		for(i=0; i<256; i++){
			fread(&bimp->colortable[i][0], 3, 1, fp);		// Reading 3 bytes into color table		
			fread(&dummy, 1, 1, fp);						// Discarding the unused byte			
		}
		for(y=0; y<bimp->ysize; y++){
			for(x=0; x<bimp->xsize; x++){
				if( fread(&colornumber, 1, 1, fp)<1 ){ sprintf(message, "Can't read pixel [%d][%d].\n", x+1, y+1); write_error_message(ErrorMessageFile, ERROR, message); }		
				bimp->image[x][y][0]=bimp->colortable[colornumber][0];		
				bimp->image[x][y][1]=bimp->colortable[colornumber][1];		
				bimp->image[x][y][2]=bimp->colortable[colornumber][2];
			}
			for(i=0; i<xzeros; i++)		if( fread(&dummy, 1, 1, fp)<1 ){ fprintf(ErrorMessageFile, "Can't read zeros after line %d.\n", y+1); write_error_message(ErrorMessageFile, ERROR, message);  }		
		}
	} else write_error_message(ErrorMessageFile, ERROR, "This is not a true color or gray scale uncompressed bitmap image. Can't read it.\n");
	fclose(fp);
	return(TRUE);
}


// Writes the specified BITMAP to file as a minimal header ascii file 

extern char out_bitmap(BITMAP *bimp, char *outfilename)
{
	unsigned char gray, dummy;
	int x, y, i, xzeros;
	FILE *fp;									// Pointer to the output file 

	if( !strcmp(outfilename, "") ) sprintf(outfilename, "%s%s~.%s", bimp->path, bimp->filename, imgext);	// Using original file name with ~ if no special name is given 
//	else { splitbmpname(outfilename, bimp); sprintf(outfilename, "%s%s.%s", bimp->path, bimp->filename, imgext); } // Adopting the given name as the new file name 

	if( (fp=fopen(outfilename, "wb"))==NULL ) { fprintf(ErrorMessageFile, "The image '%s' can't be created.\n", outfilename); exit(0); }
	if(      bimp->bitcount==24 ){ xzeros =   (bimp->xsize)%4;	                         }  // Color image:      Finding number of zeros buffered at end of each row. All rows contain 4n bytes
	else if( bimp->bitcount== 8 ){ xzeros = 4-(bimp->xsize)%4; if( xzeros==4 ) xzeros=0; }	// Grey scale image: Different number of buffered zeros
	
	dummy = 0;
	bimp->imagesize = (bimp->xsize+xzeros)*bimp->ysize;		// Updating total size of image array 
	bimp->compression=0;									// Restricting to uncompressed images
	bimp->b='B';											// Making sure bitmap gets correct signature before being written to file, for instance if read as a tiff file with signature II 
	bimp->m='M';

	fwrite(&bimp->b, 1, 1, fp);
	fwrite(&bimp->m, 1, 1, fp);
	fwrite(&bimp->filesize,       4, 1, fp);
	fwrite(&bimp->reserved,       4, 1, fp);
	fwrite(&bimp->dataoffset,     4, 1, fp);
	fwrite(&bimp->infoheadersize, 4, 1, fp);
	fwrite(&bimp->xsize,		  4, 1, fp);
	fwrite(&bimp->ysize,		  4, 1, fp);
	fwrite(&bimp->planes,		  2, 1, fp);
	fwrite(&bimp->bitcount,	 	  2, 1, fp);
	fwrite(&bimp->compression,	  4, 1, fp);
	fwrite(&bimp->imagesize,      4, 1, fp);
	fwrite(&bimp->xres,		      4, 1, fp);
	fwrite(&bimp->yres,		      4, 1, fp);
	fwrite(&bimp->colorsused,     4, 1, fp);
	fwrite(&bimp->colorsimport,   4, 1, fp);

	if( bimp->bitcount==24 ){			// It's a color image 
//		printf("Writing %d by %d pixels in RGB image to file %s.%s\n", bimp->xsize, bimp->ysize, bimp->filename, bimp->ext);
		for(y=0; y<bimp->ysize; y++){
			for(x=0; x<bimp->xsize; x++){
				fwrite(&bimp->image[x][y][0], 1, 1, fp);	// Blue  
				fwrite(&bimp->image[x][y][1], 1, 1, fp);	// Green 
				fwrite(&bimp->image[x][y][2], 1, 1, fp);	// Red   
			}
			for(i=0; i<xzeros; i++) fwrite(&dummy, 1, 1, fp);
		}
	} else if( bimp->bitcount==8 ){		// It's a gray scale image 
//		printf("Writing %d by %d pixels in gray scale image to file %s.%s\n", bimp->xsize, bimp->ysize, bimp->filename, bimp->ext);
		for(i=0; i<256; i++){			// Filling color table with a linear grayscale 
			gray = (unsigned char)i;
			fwrite(&gray,  1, 1, fp); 
			fwrite(&gray,  1, 1, fp); 
			fwrite(&gray,  1, 1, fp);
			fwrite(&dummy, 1, 1, fp); 
		}
		for(y=0; y<bimp->ysize; y++){
			for(x=0; x<bimp->xsize; x++){
				gray = (unsigned char)((bimp->image[x][y][0]+bimp->image[x][y][1]+bimp->image[x][y][2])/3.0);
				fwrite(&gray, 1, 1, fp);	// Average value of BGR channels  
			}
			for(i=0; i<xzeros; i++) fwrite(&dummy, 1, 1, fp);
		}
	}
	fclose(fp);
	return(TRUE);
}


// Copies the content of frombimp to tobimp 

extern void copy_bitmap(BITMAP *frombimp, BITMAP *tobimp)
{
	int x, y, i;
	
	strcpy(tobimp->path, frombimp->path);	
	strcpy(tobimp->filename, frombimp->filename);	
	strcpy(tobimp->ext, frombimp->ext);	
	tobimp->b=frombimp->b;
	tobimp->m=frombimp->m;
	tobimp->filesize=frombimp->filesize;
	tobimp->reserved=frombimp->reserved;
	tobimp->dataoffset=frombimp->dataoffset;
	tobimp->infoheadersize=frombimp->infoheadersize;
	tobimp->xsize=frombimp->xsize;
	tobimp->ysize=frombimp->ysize;
	tobimp->xsize4=frombimp->xsize4;
	tobimp->ysize4=frombimp->ysize4;
	tobimp->planes=frombimp->planes;
	tobimp->bitcount=frombimp->bitcount;
	tobimp->compression=frombimp->compression;
	tobimp->imagesize=frombimp->imagesize;
	tobimp->xres=frombimp->xres;
	tobimp->yres=frombimp->yres;
	tobimp->colorsused=frombimp->colorsused;
	tobimp->colorsimport=frombimp->colorsimport;

	for(y=0; y<frombimp->ysize; y++)
		for(x=0; x<frombimp->xsize; x++)
			for(i=0; i<3; i++)
				tobimp->image[x][y][i]=frombimp->image[x][y][i];
}


// Converts a gray scale bitmap to a color bitmap by simply changing the header data to those of a color image 

extern void convert2color(BITMAP *bimp)
{
	bimp->bitcount=(unsigned char)(24);
	bimp->imagesize=bimp->xsize*bimp->xsize*3;
	bimp->filesize =bimp->xsize*bimp->xsize*3+bimp->dataoffset;
	bimp->compression=0;
	bimp->dataoffset=54;
	bimp->colorsused=0;
	bimp->colorsimport=0;
}


// Crops the bitmap by the number of pixels given by the integers top, bottom, left and right 

extern char	crop_bitmap(BITMAP *bimp, int top, int bottom, int left, int right)
{
	int x, y, i;

	if( top+bottom>=bimp->ysize || left+right>=bimp->xsize ){ printf("\nError in crop_bitmap: The crop margins (top=%d bottom=%d left=%d right=%d) are larger than image (%d,%d).\n", top, bottom, left, right, bimp->xsize, bimp->ysize); return(FALSE); }

	bimp->ysize -=top;					// Cropping the top    by decreasing the image height from the top   
	bimp->xsize -=right;				// Cropping right edge by decreasing the image width  from the right 

	for(y=bottom; y<bimp->ysize; y++)	// Cropping from bottom by copying pixels down and setting new ysize 
		for(x=0; x<bimp->xsize; x++)
			for(i=0; i<3; i++) 
				bimp->image[x][y-bottom][i]=bimp->image[x][y][i];
	bimp->ysize -= bottom;

	for(x=left; x<bimp->xsize; x++)		// Cropping from left by copying pixels left and setting new xsize 
		for(y=0; y<bimp->ysize; y++)
			for(i=0; i<3; i++) 
				bimp->image[x-left][y][i]=bimp->image[x][y][i];
	bimp->xsize -= left;

	bimp->imagesize=bimp->xsize*bimp->ysize;
	bimp->filesize =bimp->imagesize+bimp->dataoffset;

	return(TRUE);
}


// Pastes the frombitmap starting from position (startx,starty) such that white is fully transparent and black is fully opaque 

extern char paste_grayscale_with_white_transparent(BITMAP *frombimp, BITMAP *tobimp, unsigned char paste_color[3], int startx, int starty)
{
	char status=TRUE;
	int x, y, i;
	int endx, endy;
	int fxs, fys;
	float transparency;

	fxs=frombimp->xsize;	// Shorthand for use below 
	fys=frombimp->ysize;
	endx = startx + fxs;
	endy = starty + fys;
	if( endx>tobimp->xsize ){ endx=tobimp->xsize; status=FALSE; }									// Making sure to never write pixels beyond the tobimp image 
	if( endy>tobimp->ysize ){ endy=tobimp->ysize; status=FALSE; }
	if( startx>tobimp->xsize || starty>tobimp->ysize ) return(FALSE);								// Making sure to never begin image beyond the tobimp image borders

	for(y=starty; y<endy; y++)
		for(x=startx; x<endx; x++)
			if( x>0 && y>0 && x-startx>0 && y-starty>0 ) {
				transparency = (float)(1.0*frombimp->image[x-startx][y-starty][0]/255); 
				for(i=0; i<3; i++) tobimp->image[x][y][i]=(unsigned char)( (1.0-transparency)*paste_color[i] + transparency*tobimp->image[x][y][i]);
			}

	return(status);
}


// Translates a color from the text string to blue-green-red values in color. Returns TRUE if the color is recognized, otherwise color is set to black 

char translate_color(char *text, unsigned char color[3])
{
	int i, l;

	l = strlen(text);
	for( i=0; i<l; i++ ) {
		if( text[i]<'a' )  text[i] += 32;						// Converting upper case letters to lower case...
		if( text[i]<'a' || text[i]>'z' ) return(FALSE);			// ...then making sure the string contains only lower case letters
	}
		 //                                      blue          green          red	
		 if( !strcmp(text, "blue"       ) ){ color[0]=255; color[1]=  0; color[2]=  0; }
	else if( !strcmp(text, "darkblue"   ) ){ color[0]=127; color[1]=  0; color[2]=  0; }
	else if( !strcmp(text, "green"      ) ){ color[0]=  0; color[1]=255; color[2]=  0; }
	else if( !strcmp(text, "darkgreen"  ) ){ color[0]=  0; color[1]=127; color[2]=  0; }
	else if( !strcmp(text, "red"        ) ){ color[0]=  0; color[1]=  0; color[2]=255; }
	else if( !strcmp(text, "darkred"    ) ){ color[0]=  0; color[1]=  0; color[2]=127; }
	else if( !strcmp(text, "cyan"       ) ){ color[0]=255; color[1]=255; color[2]=  0; }
	else if( !strcmp(text, "magenta"    ) ){ color[0]=255; color[1]=  0; color[2]=255; }
	else if( !strcmp(text, "purple"     ) ){ color[0]=127; color[1]=  0; color[2]=127; }
	else if( !strcmp(text, "yellow"     ) ){ color[0]=  0; color[1]=255; color[2]=255; }
	else if( !strcmp(text, "darkyellow" ) ){ color[0]=  0; color[1]=233; color[2]=233; }
	else if( !strcmp(text, "orange"     ) ){ color[0]=  0; color[1]=127; color[2]=255; }
	else if( !strcmp(text, "brown"      ) ){ color[0]=  0; color[1]= 63; color[2]=127; }
	else if( !strcmp(text, "black"      ) ){ color[0]=  0; color[1]=  0; color[2]=  0; }
	else if( !strcmp(text, "white"      ) ){ color[0]=255; color[1]=255; color[2]=255; }
	else if( !strcmp(text, "pink"	    ) ){ color[0]=127; color[1]=127; color[2]=255; }
	else   {								 color[0]=0;   color[1]=0;   color[2]=0;     
	//	fprintf(ErrorMessageFile, "Valid color names are: blue, green, red, darkblue, darkgreen, darkred, cyan, magenta, yellow, brown, purple, orange, black, white, pink.\n"); 
		return(FALSE); // Since not exiting here, this error message will be ignored. Better that the program doesn't stop if a non-recognized color name is given. 
	}
	return(TRUE);
}


// Splits the given textline into its path, file name and extension, returns TRUE if successfull 

char splitbmpname(char *textline, BITMAP *bimp)
{
	int i, j, l;

	l = (int)strlen(textline);
	i=l; j=0;
	while( --i>0 && textline[i]!='.' );							// Finding the first '.' from the right 
	if( i==0 ) return(FALSE);
	while( j<l-i ) bimp->ext[j]=textline[i+1+j++]; 				// Copying the letters after the last '.' to ext[] 
	bimp->ext[j]='\0';
	 textline[i]='\0';
	if(!strcmp("BMP", bimp->ext)) sprintf(bimp->ext, "bmp");	// Reducing to lower case if extension is BMP 

	l = (int)strlen(textline);									// Length minus extension 
	i=l; j=0;
	while( --i>0 && textline[i]!='\\' && textline[i]!='/');		// Finding the first '\' or '/' from the right 
	if( i>0 ){ 
		while( j<=i ) bimp->path[j]=textline[j++];				// Copying the letters from beginning until the last '\' or '/' to path[] 
		i++;													// Found a '\' or '/' 
	}
	bimp->path[j]='\0';

	while( j<=l ) bimp->filename[j-i]=textline[j++];			// Copying the rest to filename[] 
	if( strlen(bimp->filename)==0 ) return(FALSE);
	else return(TRUE);
}


// Returns the larger of the two ints				

extern int maxi(int a, int b)
{
	if( a > b ) return a;
	else return b;
}


// Returns the smaller of the two ints				

extern int mini(int a, int b)
{
	if( a < b ) return a;
	else return b;
}


// Returns the absolute value of a. Needed since overloading gave warning of truncation from double to int.				

double absd(double a)
{
	if( a < 0 ) return(-a);
	else        return( a);
}



