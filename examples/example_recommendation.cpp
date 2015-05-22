// basic file operations
#include <iostream>
#include <fstream>
#include <algorithm>    // std::max, min
#include <dai/alldai.h>  // Include main libDAI header file
#include <dai/utils/timer.h>

using namespace dai;
using namespace std;


Factor createFactorRecommendation(const Var &n1, const Var &n2, Real alpha) {
    VarSet var_set = VarSet(n1, n2);
    Factor fac(var_set);
    map<Var, size_t> state;
    for (size_t i = 0; i < n1.states(); ++i) {
        state[n1] = i;
        for (size_t j = 0; j < n2.states(); ++j) {
            state[n2] = j;
            size_t index = calcLinearState(var_set, state);

            if (i == j) {
                fac.set(index, 0.5 + alpha);
            } else {
                fac.set(index, 0.5 - alpha);
            }
        }
    }
    return fac;
}

FactorGraph example2fg() {
    vector<Var> vars;
    vector<Factor> factors;

    size_t N = 5;
    size_t M = 5;
    Real alpha = 0.1;
    // Reserve memory for the variables
    vars.reserve(N + M);

    // Create a binary variable for each movie/person
    for (size_t i = 0; i < N + M; i++)
        vars.push_back(Var(i, 2));

    factors.push_back(createFactorRecommendation(vars[0], vars[N], alpha));
    factors.push_back(createFactorRecommendation(vars[0], vars[N + 1], alpha));
    factors.push_back(createFactorRecommendation(vars[1], vars[N], alpha));
    factors.push_back(createFactorRecommendation(vars[1], vars[N + 1], alpha));
    factors.push_back(createFactorRecommendation(vars[1], vars[N + 2], alpha));
    factors.push_back(createFactorRecommendation(vars[2], vars[N + 1], alpha));
    factors.push_back(createFactorRecommendation(vars[2], vars[N + 3], alpha));
    factors.push_back(createFactorRecommendation(vars[3], vars[N + 2], alpha));
    factors.push_back(createFactorRecommendation(vars[3], vars[N + 3], alpha));
    factors.push_back(createFactorRecommendation(vars[4], vars[N + 3], alpha));
    factors.push_back(createFactorRecommendation(vars[4], vars[N + 4], alpha));
    Factor fac1(vars[5]);
    fac1.set(0, 0.9);
    fac1.set(1, 0.1);
    Factor fac2(vars[6]);
    fac2.set(0, 0.9);
    fac2.set(1, 0.1);
    factors.push_back(fac1);
    factors.push_back(fac2);

    // Create the factor graph out of the variables and factors
    cout << "Creating the factor graph..." << endl;
    FactorGraph fg(factors.begin(), factors.end(), vars.begin(), vars.end(), factors.size(), vars.size());
    string fileName = "factorGraphExampleOut.dot";
    fg.WriteToFile(fileName.c_str(), 0.01);
    return fg;
}


pair<size_t, Real> doInference(FactorGraph &fg,
                                 const string &options,
                                 size_t maxIter,
                                 Real tol,
                                 vector<Real> &m,
                                 size_t specialFactors) {
    BP bp(fg, options);

    // Choose here whether message recording should be enabled or not.
    bp.recordSentMessages = false;
    // Set properties explicitly.
    PropertySet props;
    props.set("damping", (Real)   0.0f);
    props.set("maxiter", (size_t) maxIter);
    props.set("maxtime", (Real)   INFINITY);
    props.set("tol",     (Real)   tol);
    props.set("verbose", (size_t) 0);
    props.set("specialFactors", (size_t) specialFactors);

    bp.setProperties(props);

    cout << "Properties: " << bp.printProperties() << endl;

    // Initialize inference algorithm
    //cout << "Initializing inference algorithm..." << endl;
    bp.init();

    // Initialize vector for storing the recommendations
    m = vector<Real>(fg.nrVars(), 0.0);

    // maxDiff stores the current convergence level
    Real maxDiff = 1.0;

    // Iterate while maximum number of iterations has not been
    // reached and requested convergence level has not been reached
    cout << "Running inference..." << endl;
    size_t iter;
    for (iter = 0; iter < maxIter && maxDiff > tol; iter++) {
        // Set recommendations to beliefs
        for (size_t i = 0; i < fg.nrVars(); i++)
            m[i] = bp.beliefV(i)[0];

        // Perform the requested inference algorithm for only one step
        bp.setMaxIter(iter + 1);
        maxDiff = bp.run();

        // Output progress
        //cout << "  Iterations = " << iter << ", maxDiff = " << maxDiff << endl;
    }
    
    cout << "Inference done!" << endl;
    cout << "Messages processed: " << bp.messageCount << endl;

    if (bp.recordSentMessages)
    {
        const string filename("message.txt");
        ofstream file;
        file.open(filename);
        cout << "Dumping messages to " << filename << "..." << endl;
        if (file.is_open())
            for (const auto& m : bp.getSentMessages())
                file << m.first << " <-- " << m.second << std::endl;
        file.close();
    }


    // Clean up inference algorithm

    // Return num of iterations and reached convergence level
    return make_pair(++iter, maxDiff);
}

// vector of users, containing a vector of ratings. Each rating consists of a movie id (first) and the rating (second)
vector<vector<pair<int, int> > > extract_ratings(string file_name) {
    ifstream fin;
    fin.open(file_name.c_str(), ifstream::in);
    if (!fin.is_open())
    {
        cerr << "Could not open file: " << file_name << endl;
        exit(1);
    }
    vector<vector<pair<int, int> > > ratings;
    int num_entries;
    fin >> num_entries;

    for (int i = 0; i < num_entries; ++i) {
        size_t user;
        int movie, rating;
        long long time;
        fin >> user >> movie >> rating >> time;
        user--;
        //cout << user << " " << movie << " " << rating << endl;
        while (user >= ratings.size()) {
            ratings.push_back(vector<pair<int, int> >());
        }
        ratings[user].push_back(make_pair(movie, rating));
    }
    return ratings;
}


FactorGraph data2fg(const vector<vector<pair<int, int> > > &votings, int user, size_t& specialFactors) {
    // We will create a variable for every potential user/movie. We know the number of users, let us estimate the number of movies.
    int num_users = votings.size();
    int num_movies = 0;
    for (size_t i = 0; i < votings.size(); ++i) {
        if (votings[i].size() > 0) {
            num_movies = max(num_movies, votings[i][votings[i].size() - 1].first);
        }
    }
    // We add one to avoid the case where movies start counting at 1, leading to one additional movie.
    num_movies++;

    Real alpha = 0.0001;
    int threshold = 4;
    vector<Var> vars;
    vector<Factor> factors;

    // Reserve memory for the variables
    cout << "Estimated num_users/num_movies: " << num_users << "/" << num_movies << endl;
    vars.reserve(num_users + num_movies);
    // Create a binary variable for each movie/person
    for (size_t i = 0; i < (size_t) (num_users + num_movies); i++)
        vars.push_back(Var(i, 2));

    for (size_t i = 0; i < votings.size(); ++i) {
        for (size_t j = 0; j < votings[i].size(); ++j) {
            if (votings[i][j].second >= threshold) {
                factors.push_back(createFactorRecommendation(vars[i], vars[num_users + votings[i][j].first], alpha));
            }
        }
    }
    specialFactors = factors.size();

    cout << "Factors created. Dealing with the user now..." << endl;
    // calculate some metrics for the user.
    int normalization_factor_p = 4;
    Real sum = 0;
    Real sq_sum = 0;
    for (size_t i = 0; i < votings[user].size(); ++i) {
        sum += votings[user][i].second;
    }
    Real mean = sum / votings[user].size();
    for (size_t i = 0; i < votings[user].size(); ++i) {
        sq_sum += pow(votings[user][i].second - mean, 2);
    }
    Real stdev = std::sqrt(sq_sum / votings[user].size());

    for (size_t i = 0; i < votings[user].size(); ++i) {
        Factor fac(vars[num_users + votings[user][i].first]);
        Real like = max(0.1, min(0.9, 0.5 + (votings[user][i].second - mean) / (stdev * normalization_factor_p)));
        //cout << votings[user][i].second << " to " << like << endl;
        fac.set(0, like);
        fac.set(1, 1-like);
        factors.push_back(fac);
    }

    // Create the factor graph out of the variables and factors
    cout << "Creating the factor graph..." << endl;
    return FactorGraph(factors.begin(), factors.end(), vars.begin(), vars.end(), factors.size(), vars.size());
}

pair<Real, Real> getPrecisionAndRecall(const vector<vector<pair<int, int> > > &test_data,
                                           const vector<pair<Real, int> > &ratings, int user, size_t N) {
    // get the top predicted elements.
    vector<int> predicted;
    for (size_t i = 0; i < std::min(N, ratings.size()); ++i) {
        predicted.push_back(ratings[i].second);
    }

    vector<int> wanted;
    for (size_t i = 0; i < test_data[user].size(); ++i) {
        if (test_data[user][i].second >= 5) {
            wanted.push_back(test_data[user][i].first);
        }
    }

    // compute the intersection.
    vector<int> v(N);
    sort(predicted.begin(), predicted.end());
    sort(wanted.begin(), wanted.end());
    vector<int>::iterator it = set_intersection(predicted.begin(), predicted.end(), wanted.begin(), wanted.end(),
                                                v.begin());
    v.resize(it - v.begin());

    //std::cout << "The intersection has " << (v.size()) << " elements:\n";
    if (wanted.size() == 0) {
        // dummy element to prevent division by zero.
        wanted.push_back(0);
    }
    return make_pair<Real, Real>(v.size() / static_cast<Real>(N), v.size() / static_cast<Real>(wanted.size()));
}

// compares the calculated ratings to the reference file
void verify_results(const string dataset, const Real delta, vector<pair<Real, int>>& ratings){
    ifstream fin;
    string file_name = dataset + ".reference";
    fin.open(file_name.c_str(), ifstream::in);
    if(!fin.is_open())
    {
        cerr << "Test FAILED! Reference file " << dataset << ".reference cannot be opened." << endl;
        exit(1);
    }
    Real ref_d;
    Real ref_i;
    size_t line = 1;
    for(vector<pair<Real, int> >::iterator it=ratings.begin();it!=ratings.end();it++,line++) {
        fin >> ref_d;
        fin >> ref_i;
        if(fabs(ref_d - it->first) > delta || ref_i != it->second)
        {
            cerr << "TEST FAILED! Difference at line " << line << " found!" << endl;
            cerr << "Expected: " << ref_d << " / " << ref_i << endl;
            cerr << "Actual: " << it->first << " / " << it->second << endl;
            exit(1);
        }
    }
    cerr << "TEST SUCCESSFULL!" << endl;
    fin.close();
}

/// Main program
int main(int argc, char **argv) {
    const size_t maxiter =  100;
    const Real tol = 1e-15;
    const string options("[updates=SEQMAX,maxiter=100,tol=1e-15,logdomain=0]");

    // default cpufreq to 2.6Ghz
    const long cpu_freq =26e+8;
    const bool output_ratings = false;
    const string dataset ="u1";
    const bool run_tests = false;
    const Real delta = 1e-12;
    const int num_measurements =  1;

    cout << "reading " << dataset << ".base now..." << endl;
    vector<vector<pair<int, int> > > input_data = extract_ratings(dataset + ".base");
    vector<vector<pair<int, int> > > test_data = extract_ratings(dataset + ".test");
    const int N = 1;
    Real p10 = 0;
    Real p20 = 0;
    Real r10 = 0;
    Real r20 = 0;
    long measured_cycles = 0;
    Timer timer;
    vector<long> measurements;
    for(int run=0;run<num_measurements;run++) {
        p10 = 0;
        p20 = 0;
        r10 = 0;
        r20 = 0;
        measured_cycles = 0;
        for (int user=0; user<N; ++user) {
            cout << "building factor graph for user " << user+1 << " out of " << N << endl;
            //FactorGraph fg  = example2fg();
            size_t specialFactors = -1;
            FactorGraph fg = data2fg(input_data, user, specialFactors);

            vector<Real> m; // Stores the final recommendations
            cout << "Inference algorithm: BP. Options: " << options << endl;
            cout << "Solving the inference problem...please be patient!" << endl;
            cout << "Note: There's no output during the inference. You may have to wait a bit..." << endl;

            timer.tic();
            pair<size_t, Real> result = doInference(fg, options, maxiter, tol, m, specialFactors);
            measured_cycles += timer.toc();

            cout << "Iterations = " << result.first << ", maxDiff = " << result.second << endl;

            // Debug output with all the values:
            // for (size_t i = 0; i < m.size(); ++i) {
            //    cout << i << " with " << m[i] << endl;
            // }

            vector<pair<Real, int> > ratings;
            for (size_t i = input_data.size(); i < m.size(); ++i) {
                // push back the negative so we can use the standard sorting.
                ratings.push_back(make_pair<Real, int>(-m[i], i - input_data.size() + 1));
            }
            cout << "RatingsSize" << ratings.size() << endl;
            if(output_ratings && run == num_measurements -1) {
                // output the calculated ratings to STDERR so that they can be stored and reused for regression tests
                // you can create a reference file the following way:
                //      ./example_recommendation > output.txt 2> ratings.txt
                cerr.precision(15);
                for(vector<pair<Real, int> >::iterator it=ratings.begin();it!=ratings.end();it++) {
                    cerr << it->first << " " << it->second << endl;
                }
            }
            if(run_tests)
                verify_results(dataset, delta, ratings);

            sort(ratings.begin(), ratings.end());
            pair<Real, Real> pr10 = getPrecisionAndRecall(test_data, ratings, user, 10);
            pair<Real, Real> pr20 = getPrecisionAndRecall(test_data, ratings, user, 20);

            p10 +=  pr10.first;
            p20 +=  pr20.first;
            r10 +=  pr10.second;
            r20 +=  pr20.second;
            cout << "Precision (N=10): " << pr10.first << endl;
            cout << "Precision (N=20): " << pr20.first << endl;
            cout << "Recall (N=10): " << pr10.second << endl;
            cout << "Recall (N=20): " << pr20.second << endl;
            measurements.push_back(measured_cycles);
        }
    }
    p10 = p10 / static_cast<Real>(N);
    p20 = p20 / static_cast<Real>(N);
    r10 = r20 / static_cast<Real>(N);
    r20 = r20 / static_cast<Real>(N);
    sort(measurements.begin(), measurements.end());
    measured_cycles = measurements[num_measurements / 2];
    cout << "Final estimated:" << endl;
    cout << "Precision (N=10): " << p10 << endl;
    cout << "Precision (N=20): " << p20 << endl;
    cout << "Recall (N=10): " << r10 << endl;
    cout << "Recall (N=20): " << r20 << endl;
    cout << "Ran " << num_measurements << " times.\nMedian:" << endl;
    cout << "Measured cycles: " << measured_cycles << endl;
    cout << "Runtime: " << ((double) measured_cycles) / cpu_freq << " seconds" << endl;
    return 0;
}
