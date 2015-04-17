#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>

using namespace std;

void extract_movies(map<int, int>& found, int user_id, size_t max_movies, string file_name) {
    ifstream fin;
    fin.open(file_name.c_str(), ifstream::in);
    vector<vector<pair<int, int> > > ratings;
    int num_entries;
    fin >> num_entries;
    for (int i = 0; i < num_entries; ++i) {
        int user, movie, rating;
        long long time;
        fin >> user >> movie >> rating >> time;
        if (user == user_id) {
        	if (found.size() < max_movies) {
	        	found[movie] = found.size();
	        }
        }
    }
}

vector<tuple<int, int, int, long long> > remove_unused(map<int, int>& found, string file_name1) {
    ifstream fin;
    fin.open(file_name1.c_str(), ifstream::in);
    vector<tuple<int, int, int, long long> > ratings;
    int num_entries;
    fin >> num_entries;
    for (int i = 0; i < num_entries; ++i) {
        int user, movie, rating;
        long long timeT;
        fin >> user >> movie >> rating >> timeT;
        if (found.find(movie) != found.end() )
		{
		    ratings.push_back(make_tuple(user, found[movie], rating, timeT));
		}
    }
    return ratings;
}


void write_file(vector<tuple<int, int, int, long long> >& ratings, string file_name1) {
    ofstream fout;
    fout.open(file_name1.c_str(), ofstream::out);
    fout << ratings.size() << endl;
    for (size_t i = 0; i < ratings.size(); ++i) {
        fout << get<0>(ratings[i]) << "\t"<< get<1>(ratings[i]) << "\t"<< get<2>(ratings[i]) << "\t"<< get<3>(ratings[i]) << endl;
    }
}



int main(void) {
	const int user_id = 1;
	const size_t max_movies = 50;
	map<int, int> found;
	extract_movies(found, user_id, max_movies, "u1.base");
	extract_movies(found, user_id, max_movies, "u1.test");

	vector<tuple<int, int, int, long long> > ratings1 =remove_unused(found, "u1.base");
	vector<tuple<int, int, int, long long> > ratings2 =remove_unused(found, "u1.test");
	write_file(ratings1, "uV2New"+to_string(user_id)+".base");
	write_file(ratings2, "uV2New"+to_string(user_id)+".test");

	return 0;
}
