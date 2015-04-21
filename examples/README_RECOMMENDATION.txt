example_recommendation supports the following additional commandline arguments:

-cpu_freq 
    Can be used to specify the frequency (for correctly calculating the runtime in seconds)

-printRatings
    Prints the calculated ratings to STDERR (!)

-dataset
    Specify which file in ml-100k/ to use as input. No file-extension needed. e.g "-dataset u1"

-test
    Compare the calculated ratings with the reference file. We fail with exit(1) if there's a difference

-testDelta
   A test fails if fabs(rating - reference_rating) > testDelta.  



EXAMPLES
# create a new reference file
./example_recommendation -dataset u1 -printRatings 2> u1.referece

# run tests and output cycles and seconds
./example_recommendation -dataset u1 -test | tail -n 2

