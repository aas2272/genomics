README

There are two main files used in this project.  The first is run.py, it contains the main method to run the analysis.  The input is a set of tissue type samples with underscores between words.  The last input is a number representing which basis element should be displayed in the output charts.  For example, the like run when we did it was:

python run.py Blood_Whole_Blood Lung_Lung Thyroid_Thyroid Ovary_Ovary Nerve_Nerve_-_Tibial 3

The first tissue type is where the miniumum variance is drawn from.

The main method reads the first tissue and either gets the already parsed and saved chart, if it exists, or creates it from the data.gct file what was taken from the gtex portal (this data needs to be downloaded directly from the gtex website).  It then cuts down to only the highly expressed genes (70% above .1) and then gets all the variances.  We then get a bunch of different variances (smallest largest percentile, median, mean, etc...)

We then parse the data for just the columns (samples) that correspond to the 5 tissue types inputted, and then re-cut the rows (genes to the prevalent genes).

Using one of the variances we create a matrix of the same size as the cut down data matrix with entries normally distributed (0, sqrt(var)).  This variance is also used as the cutoff for the the difference between the random normal matrix and our factorized matrix.

We then factorize the data matrix to matricies rank 1,2,3,... until the the corresponding product is within the cutoff of the original data matrix.  This is done using the find lim function with the get diff function.  We rank we found was 4.

Lastly we create a chart with the data.  We go through each input tissue, and then through each of the basis elements and then through each of the samples of that tissue type, continuously summing up the total (to get the average), and also keeping track of the values to make the histogram.  The equation for the sume is given by the function called fraction.  We save this output and create histograms out of the arrays.

sample_parser.py is the other main file.  It keeps track and manages the data and can stand alone.

It is initialized with the data to parse, whether this tissue has been parsed already (1 if its new 2 if it was already created), the map of sample names to tissue types, the tissue types to parse.

It has the following methods:
get_samples: returns the types of tissues

is_valid: determines is the desired tissues exist in the data

preprocess: stand alone method to get the data out of the file

save_file: safes the data to a file

get_data: reads in the data to be used with this class

parse_genes: reads the mapping and saves them in a dictionary

get_relevant_tissues: cuts the data down to only the samples of the tissues desired

get_relevant_genes: cuts the data down to only genes that are heavily expressed

get_var: returns a list of all variances

Additionally, the full data needs to be stored as data.gct and the tissue map is tissues.txt both of which are found on the gtex portal website and there needs to be a data folder where all the intermediate data goes.

requirements.txt has all of the requirements needed to run the code

