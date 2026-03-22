#include <bits/stdc++.h>
#include <span>
#include <filesystem>
#include "graph.h"
using namespace std;
namespace fs = std::filesystem;

int main() {
    string directory = "../data/graphs/infoarena/";
    
    // Iterate through all grader_test*.in files
    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.is_regular_file()) {
            string filename = entry.path().filename().string();
            
            // Check if file matches pattern grader_test#.in
            if (filename.find("grader_test") == 0 && filename.find(".in") != string::npos) {
                string input_path = entry.path().string();
                
                // Extract test number from filename
                size_t start = filename.find("grader_test") + 11;
                size_t end = filename.find(".in");
                string test_num = filename.substr(start, end - start);
                
                string output_path = directory + "test" + test_num + ".txt";
                
                // Open input file
                ifstream fin(input_path);
                if (!fin) {
                    cerr << "Could not open input file: " << input_path << endl;
                    continue;
                }
                
                NodeID n;
                EdgeID m;
                fin >> n >> m;
                
                vector<tuple<NodeID, NodeID, Distance>> edges;
                for (EdgeID i = 0; i < m; i++) {
                    NodeID x, y;
                    Distance c;
                    fin >> x >> y >> c;
                    x--, y--;
                    
                    // Create edge from x to y with weight c
                    edges.emplace_back(x, y, c);
                }
                
                fin.close();
                
                // Create graph from edges
                Graph g(n, edges);
                
                // Print graph to output file
                ofstream fout(output_path);
                if (!fout) {
                    cerr << "Could not open output file: " << output_path << endl;
                    continue;
                }
                
                g.format_print(fout, false);
                fout.close();
                
                cout << "Converted " << input_path << " -> " << output_path << endl;
            }
        }
    }
    
    return 0;
}