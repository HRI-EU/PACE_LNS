/*
Copyright (c) 2024
Honda Research Institute Europe GmbH
Carl-Legien-Str. 30
63073 Offenbach/Main
Germany

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <map>
#include <math.h>
#include <numeric>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <vector>

using namespace std;
    
struct Node {
  int idx;
  double order_val;
  bool operator==(const Node &rhs) const { return idx == rhs.idx; }
  bool operator<(const Node &rhs) const { return order_val < rhs.order_val; }
  bool operator>(const Node &rhs) const { return order_val > rhs.order_val; }
};

class Problem {
public:
  int n_left, n_right, n_edges;
  std::vector<Node> r_nodes;
  double exact_opt = -1.0;
  std::vector<std::vector<bool>> adj_matrix;
  std::vector<std::vector<int>> crossing_matrix;
  std::vector<std::vector<int>> crossing_matrix_transposed;
  std::vector<std::vector<int>> graph;

  std::vector<int> crossings_first_position;
  void read_in(std::string filename,
               int max_nodes_crossing_computation = 15000) {
    graph.clear();
    exact_opt = -1.0;
    r_nodes = {};
    adj_matrix = {};
    crossing_matrix = {};
    crossing_matrix_transposed = {};
    std::ifstream in(filename);
    if (in.is_open()) {
      string trash_s;
      int n_l, n_r, trash;
      in >> trash_s >> trash_s >> n_left >> n_right >> n_edges;
      adj_matrix = std::vector<std::vector<bool>>(
          n_left, std::vector<bool>(n_right, false));

      for (int n = 1; n <= n_right; n++) {
        graph.push_back({});

        Node node;
        node.idx = n_left + n;
        r_nodes.push_back(node);
      }
      for (int e = 0; e < n_edges; e++) {
        in >> n_l >> n_r;

        graph[n_r - n_left - 1].push_back(n_l);

        adj_matrix[n_l - 1][n_r - (n_left + 1)] = true;
      }
      // we need the neighbors of the right-side nodes sorted, which is not the
      // case in the public instances
      for (int n = 1; n <= n_right; n++) {
        std::sort(graph[n - 1].begin(), graph[n - 1].end());
      }
    } else {
      printf("Cannot open input file %s\n", filename.c_str());
      exit(1);
    }
    if (n_right <= max_nodes_crossing_computation) {
      compute_crossing_matrix();
    }
  }

  void read_from_cin(int max_nodes_crossing_computation = 15000) {
    graph.clear();
    exact_opt = -1.0;
    r_nodes = {};
    adj_matrix = {};
    crossing_matrix = {};
    crossing_matrix_transposed = {};

    string trash_s;
    int n_l, n_r, trash;
    ios_base::sync_with_stdio(false);
    cin.tie();
    cin >> trash_s >> trash_s >> n_left >> n_right >> n_edges;

    adj_matrix = std::vector<std::vector<bool>>(
        n_left, std::vector<bool>(n_right, false));

    for (int n = 1; n <= n_right; n++) {
      graph.push_back({});

      Node node;
      node.idx = n_left + n;
      r_nodes.push_back(node);
    }
    for (int e = 0; e < n_edges; e++) {
      cin >> n_l >> n_r;

      graph[n_r - n_left - 1].push_back(n_l);

      adj_matrix[n_l - 1][n_r - (n_left + 1)] = true;
    }
    // we need the neighbors of the right-side nodes sorted, which is not the
    // case in the public instances
    for (int n = 1; n <= n_right; n++) {
      std::sort(graph[n - 1].begin(), graph[n - 1].end());
    }

    if (n_right <= max_nodes_crossing_computation) {
      compute_crossing_matrix();
    }
  }

  Problem get_subproblem(std::vector<int> nodes) {
    Problem subproblem;
    subproblem.n_right = (int)nodes.size();
    subproblem.n_left = n_left;
    subproblem.adj_matrix = std::vector<std::vector<bool>>(
        subproblem.n_left, std::vector<bool>(subproblem.n_right, false));
    subproblem.graph = {};
    for (int n = 0; n < nodes.size(); n++) {
      subproblem.graph.push_back(graph[nodes[n] - n_left - 1]);
      for (int l = 0; l < graph[nodes[n] - n_left - 1].size(); l++) {
        subproblem.adj_matrix[graph[nodes[n] - n_left - 1][l] - 1][n] = true;
      }
      Node node;
      node.idx = n_left + n + 1;
      subproblem.r_nodes.push_back(node);
    }
    // note: we do not compute adj_matrix here, so segtree eval cannot be used
    // on subproblem

    // subproblem.compute_crossing_matrix();

    return subproblem;
  }

  void compute_crossing_matrix() {
    crossing_matrix =
        std::vector<std::vector<int>>(n_right, std::vector<int>(n_right, 0));
    crossing_matrix_transposed =
        std::vector<std::vector<int>>(n_right, std::vector<int>(n_right, 0));
    for (int i = 0; i < n_right; i++) {
      for (int j = i + 1; j < n_right; j++) {
        int n_lower = 0;
        int n_upper = 0;
        if (graph[i].size() < graph[j].size()) {
          for (int k = 0; k < graph[i].size(); k++) {
            if (n_lower < graph[j].size()) {
              // number of elements in graph[j], which are <= graph[i][k]
              n_lower = (int)(upper_bound(graph[j].begin() + n_lower,
                                          graph[j].end(), graph[i][k]) -
                              graph[j].begin());
              n_upper = (int)graph[j].size() - n_lower;
              // n_lower has to be decreased if graph[i][k] is in graph[j]
              if (n_lower > 0 && graph[j][n_lower - 1] == graph[i][k]) {
                n_lower--;
              }
            }
            crossing_matrix[i][j] += n_lower;
            crossing_matrix[j][i] += n_upper;

            crossing_matrix_transposed[j][i] += n_lower;
            crossing_matrix_transposed[i][j] += n_upper;
          }
        } else {
          for (int k = 0; k < graph[j].size(); k++) {
            if (n_lower < graph[i].size()) {
              // number of elements in graph[i], which are <= graph[j][k]
              n_lower = (int)(upper_bound(graph[i].begin() + n_lower,
                                          graph[i].end(), graph[j][k]) -
                              graph[i].begin());
              n_upper = (int)graph[i].size() - n_lower;
              // n_lower has to be decreased if graph[j][k] is in graph[i]
              if (n_lower > 0 && graph[i][n_lower - 1] == graph[j][k]) {
                n_lower--;
              }
            }
            crossing_matrix[j][i] += n_lower;
            crossing_matrix[i][j] += n_upper;

            crossing_matrix_transposed[i][j] += n_lower;
            crossing_matrix_transposed[j][i] += n_upper;
          }
        }
      }
    }

    for (int i = 0; i < n_right; i++) {
      crossings_first_position.push_back(std::accumulate(
          crossing_matrix[i].begin(), crossing_matrix[i].end(), 0));
    }
  }

  void random_init(int seed, int l_min, int l_max, int r_min, int r_max,
                   double p_min, double p_max,
                   int max_nodes_crossing_computation = 15000) {
    graph.clear();
    exact_opt = -1.0;
    r_nodes = {};
    adj_matrix = {};
    crossing_matrix = {};
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> n_left_dist(l_min, l_max);
    std::uniform_int_distribution<int> n_right_dist(r_min, r_max);
    std::uniform_real_distribution<double> prob_dist(p_min, p_max);
    std::uniform_real_distribution<double> p_dist(0.0, 1.0);

    n_left = n_left_dist(gen);
    n_right = n_right_dist(gen);
    adj_matrix = std::vector<std::vector<bool>>(
        n_left, std::vector<bool>(n_right, false));
    double prob = prob_dist(gen);
    double p;
    for (int i = n_left + 1; i < n_left + n_right + 1; i++) {
      graph[i - n_left - 1] = {};
      Node node;
      node.idx = i;
      r_nodes.push_back(node);
      for (int j = 1; j < n_left + 1; j++) {
        p = p_dist(gen);
        if (p < prob) {
          adj_matrix[j - 1][i - (n_left + 1)] = true;
          graph[i - n_left - 1].push_back(j);
        }
      }
    }
    if (n_right <= max_nodes_crossing_computation) {
      compute_crossing_matrix();
    }
  }

  void write(std::string filename) {
    ofstream outfile;
    outfile.open(filename);
    outfile << "p ocr " << n_left << " " << n_right << " " << n_edges << "\n";
    for (int l = 0; l < n_left; l++) {
      for (int r = 0; r < n_right; r++) {
        if (adj_matrix[l][r]) {
          outfile << l + 1 << " " << r + n_left + 1 << "\n";
        }
      }
    }
    outfile.close();
  }

  void write_solution(std::string filename, std::vector<int> solution) {
    ofstream outfile;
    outfile.open(filename);
    for (int i = 0; i < solution.size(); i++) {
      outfile << solution[i] << "\n";
    }
    outfile.close();
  }

  double eval_trivial(std::vector<int> &solution) {
    double crossings = 0.0;
    int n1, n2, l1, l2;
    for (int i = 0; i < solution.size(); i++) {
      n1 = solution[i];
      for (int j = 0; j < graph[n1 - n_left - 1].size(); j++) {
        l1 = graph[n1 - n_left - 1][j];
        for (int k = i + 1; k < solution.size(); k++) {
          n2 = solution[k];
          for (int l = 0; l < graph[n2 - n_left - 1].size(); l++) {
            l2 = graph[n2 - n_left - 1][l];
            if (l2 < l1) {
              crossings += 1.0;
            } else {
              // we can finish here since nodes in graph are sorted
              break;
            }
          }
        }
      }
    }
    if (exact_opt == -1.0) {
      return crossings;
    } else {
      return crossings - exact_opt;
    }
  }

  // helper function for segment tree
  long sum(std::vector<long> &segment_tree, int L, int R, int sl, int sr,
           int i) {
    if (sl >= R || sr <= L)
      return 0; // Outside of range
    if (sl >= L && sr <= R)
      return segment_tree[i]; // Inside of range
    int mid = (sl + sr) / 2;
    return sum(segment_tree, L, R, sl, mid, i * 2 + 1) +
           sum(segment_tree, L, R, mid, sr, i * 2 + 2);
  }

  // evaluation with segment tree based on PACE code and tutorial at
  // https://blog.garybricks.com/segment-tree-introduction-in-c
  double eval_segment(std::vector<int> &solution) {
    long crossings = 0;
    int size = (int)solution.size();
    int base_size = 1;
    int seg_tree_size;
    if (size <= 1) {
      if (exact_opt == -1.0) {
        return 0.0;
      } else {
        return -exact_opt;
      }
    }

    std::vector<int> right_order(solution.size());
    for (int i = 0; i < solution.size(); i++) {
      right_order[solution[i] - n_left - 1] = i;
    }
    while (base_size < size) {
      base_size *= 2;
    }
    seg_tree_size = base_size * 2 - 1;
    std::vector<long> segtree(seg_tree_size, 0);

    for (int l = 0; l < n_left; l++) {
      for (int r = 0; r < solution.size(); r++) {
        int nr = solution[r];
        if (adj_matrix[l][nr - (n_left + 1)]) {
          int i = right_order[nr - n_left - 1] + seg_tree_size / 2;
          segtree[i] += 1;
          while (i > 0) {
            i = (i - 1) / 2;
            segtree[i] = segtree[i * 2 + 1] + segtree[i * 2 + 2];
          }
          if (size > right_order[nr - n_left - 1] + 1) {
            crossings += sum(segtree, right_order[nr - n_left - 1] + 1, size, 0,
                             base_size, 0);
          }
        }
      }
    }

    if (exact_opt == -1.0) {
      return (double)crossings;
    } else {
      return (double)crossings - exact_opt;
    }
  }

  double eval(std::vector<int> &solution) {
    if (crossing_matrix.size() == 0) {
      return eval_segment(solution);
    }
    long crossings = 0;
    if ((int)solution.size() <= 1) {

      if (exact_opt == -1.0) {
        return 0.0;
      } else {
        return -exact_opt;
      }
    }

    for (int i = 0; i < solution.size(); i++) {
      for (int j = i + 1; j < solution.size(); j++) {
        crossings += crossing_matrix[solution[i] - (n_left + 1)]
                                    [solution[j] - (n_left + 1)];
      }
    }

    if (exact_opt == -1.0) {
      return (double)crossings;
    } else {
      return (double)crossings - exact_opt;
    }
  }
};

std::vector<int> random_insertion(Problem &problem, std::vector<int> &solution,
                                  std::vector<int> &new_nodes, double &obj_out,
                                  std::mt19937 &gen, double obj_start = -1.0) {
  std::vector<int> result = solution;
  std::uniform_int_distribution<int> dist;
  double obj = 0.0;
  if (obj_start >= 0.0) {
    obj = obj_start;
  }
  for (int n = 0; n < new_nodes.size(); n++) {
    dist = std::uniform_int_distribution<int>(0, result.size());
    int node = new_nodes[n];
    int pos = dist(gen);
    if (obj_start >= 0.0 &&
        problem.graph[node - problem.n_left - 1].size() > 0) {
      for (int i = 0; i < pos; i++) {
        obj +=
            (double)problem
                .crossing_matrix_transposed[node - (problem.n_left + 1)]
                                           [result[i] - (problem.n_left + 1)];
      }
      for (int i = pos; i < result.size(); i++) {
        obj +=
            (double)problem.crossing_matrix[node - (problem.n_left + 1)]
                                           [result[i] - (problem.n_left + 1)];
      }
    }
    result.insert(result.begin() + pos, node);
  }
  if (obj_start >= 0.0) {
    obj_out = obj;
  } else {
    obj_out = problem.eval(result);
  }
  return result;
}

std::vector<int> greedy_insertion(Problem &problem, std::vector<int> &solution,
                                  std::vector<int> &new_nodes, double &obj_out,
                                  double obj_start = -1.0) {
  std::vector<int> result = solution;
  double obj;
  if (obj_start >= 0.0) {
    obj = obj_start;
  } else {
    obj = problem.eval(solution);
  }

  for (int n = 0; n < new_nodes.size(); n++) {
    int node = new_nodes[n];

    if (result.size() == 0) {
      result.push_back(node);
      continue;
    }
    int size = (int)result.size() + 1;

    int best_pos;
    long best_cross;
    if (problem.crossing_matrix.size() > 0) {
      int c_before = 0;
      int c_after = 0;
      int cross;
      int cross_noise;

      // first compute #crossings if node is at beginning of sequence and
      // substract crossings with nodes currently not in solution
      c_after = problem.crossings_first_position[node - (problem.n_left + 1)];
      if (c_after > 0) {
        for (int n2 = n + 1; n2 < new_nodes.size(); n2++) {
          c_after -=
              problem.crossing_matrix[node - (problem.n_left + 1)]
                                     [new_nodes[n2] - (problem.n_left + 1)];
        }
      }

      best_cross = c_after;
      best_pos = 0;
      int idx = node - (problem.n_left + 1);
      int i = 0;
      for (int node2 : result) {
        c_after -= problem.crossing_matrix[idx][node2 - (problem.n_left + 1)];
        c_before +=
            problem
                .crossing_matrix_transposed[idx][node2 - (problem.n_left + 1)];

        cross = c_before + c_after;

        if (c_before >= best_cross) {
          break;
        }
        if (cross < best_cross) {
          best_cross = cross;
          best_pos = i + 1;
          if (c_after == 0 || best_cross == 0) {
            break;
          }
        }
        ++i;
      }

    } else {
      best_pos = 0;
      best_cross = 0;
      if (result.size() > 0) {
        std::vector<int> before((int)result.size(), 0);
        std::vector<int> after((int)result.size(), 0);
        for (int j = 0; j < result.size(); j++) {
          int node2 = result[j];
          for (int k = 0; k < problem.graph[node - problem.n_left - 1].size();
               k++) {
            // number of elements in graph[n_left+node2+1], which are >=
            // graph[n_left+node+1][k]
            int n_lower =
                (int)(upper_bound(
                          problem.graph[node2 - problem.n_left - 1].begin(),
                          problem.graph[node2 - problem.n_left - 1].end(),
                          problem.graph[node - problem.n_left - 1][k]) -
                      problem.graph[node2 - problem.n_left - 1].begin());
            int n_upper =
                (int)problem.graph[node2 - problem.n_left - 1].size() - n_lower;

            // n_lower has to be decreased if graph[n_left+i+1][k] is in
            // graph[n_left+j+1]
            if (n_lower > 0 &&
                problem.graph[node2 - problem.n_left - 1][n_lower - 1] ==
                    problem.graph[node - problem.n_left - 1][k]) {
              n_lower--;
            }
            before[j] +=
                n_lower; // number of crossings if node is placed before node2
            after[j] +=
                n_upper; // number of crossings if node is placed after node2
          }
        }
        long cross = std::accumulate(before.begin(), before.end(), 0);
        best_cross = cross;
        for (int i = 1; i < result.size() + 1; i++) {
          cross = cross - before[i - 1] + after[i - 1];
          if (cross < best_cross) {
            best_pos = i;
            best_cross = cross;
          }
        }
      }
    }
    result.insert(result.begin() + best_pos, node);
    obj += (double)best_cross;
  }

  obj_out = obj;

  return result;
}

std::vector<int> barycenter(Problem &problem) {
  std::vector<Node> nodes = problem.r_nodes;
  for (int n = 0; n < nodes.size(); n++) {
    if (problem.graph[nodes[n].idx - problem.n_left - 1].size() == 0) {
      nodes[n].order_val = 0.0;
    } else {
      nodes[n].order_val =
          accumulate(problem.graph[nodes[n].idx - problem.n_left - 1].begin(),
                     problem.graph[nodes[n].idx - problem.n_left - 1].end(),
                     0.0) /
          problem.graph[nodes[n].idx - problem.n_left - 1].size();
    }
  }
  std::sort(nodes.begin(), nodes.end());
  std::vector<int> result;
  for (int n = 0; n < nodes.size(); n++) {
    result.push_back(nodes[n].idx);
  }
  return result;
}

std::vector<int> select_block(Problem &problem, std::vector<int> &solution,
                              int k, std::mt19937 &gen, double &obj_start) {
  std::vector<int> new_nodes;

  std::uniform_int_distribution<int> dist =
      std::uniform_int_distribution<int>(0, solution.size() - 1 - k);
  int idx = dist(gen);
  for (int i = 0; i < k; i++) {
    new_nodes.push_back(solution[idx]);
    // substract crossings of the node to remove from objective value
    if (problem.graph[solution[idx] - problem.n_left - 1].size() > 0) {
      /*for(int i=0; i<idx; i++){
              //obj_start -=
      (double)problem.crossing_matrix[solution[i]-(problem.n_left+1)][solution[idx]-(problem.n_left+1)];
              obj_start -=
      (double)problem.crossing_matrix_transposed[solution[idx]-(problem.n_left+1)][solution[i]-(problem.n_left+1)];
      }

      for(int i=idx+1; i<solution.size(); i++){
              obj_start -=
      (double)problem.crossing_matrix[solution[idx]-(problem.n_left+1)][solution[i]-(problem.n_left+1)];
      }*/

      // unrolled version of the two previous loops
      int repeat = idx / 5;
      int left = idx % 5;
      int j = 0;
      while (repeat--) {
        obj_start -= (double)(problem.crossing_matrix_transposed
                                  [solution[idx] - (problem.n_left + 1)]
                                  [solution[j] - (problem.n_left + 1)] +
                              problem.crossing_matrix_transposed
                                  [solution[idx] - (problem.n_left + 1)]
                                  [solution[j + 1] - (problem.n_left + 1)] +
                              problem.crossing_matrix_transposed
                                  [solution[idx] - (problem.n_left + 1)]
                                  [solution[j + 2] - (problem.n_left + 1)] +
                              problem.crossing_matrix_transposed
                                  [solution[idx] - (problem.n_left + 1)]
                                  [solution[j + 3] - (problem.n_left + 1)] +
                              problem.crossing_matrix_transposed
                                  [solution[idx] - (problem.n_left + 1)]
                                  [solution[j + 4] - (problem.n_left + 1)]);
        j += 5;
      }
      for (int l = 0; l < left; l++) {
        obj_start -= (double)problem.crossing_matrix_transposed
                         [solution[idx] - (problem.n_left + 1)]
                         [solution[j + l] - (problem.n_left + 1)];
      }

      j = idx + 1;
      repeat = ((int)solution.size() - j) / 5;
      left = ((int)solution.size() - j) % 5;
      while (repeat--) {
        obj_start -=
            (double)(problem
                         .crossing_matrix[solution[idx] - (problem.n_left + 1)]
                                         [solution[j] - (problem.n_left + 1)] +
                     problem
                         .crossing_matrix[solution[idx] - (problem.n_left + 1)]
                                         [solution[j + 1] -
                                          (problem.n_left + 1)] +
                     problem
                         .crossing_matrix[solution[idx] - (problem.n_left + 1)]
                                         [solution[j + 2] -
                                          (problem.n_left + 1)] +
                     problem
                         .crossing_matrix[solution[idx] - (problem.n_left + 1)]
                                         [solution[j + 3] -
                                          (problem.n_left + 1)] +
                     problem
                         .crossing_matrix[solution[idx] - (problem.n_left + 1)]
                                         [solution[j + 4] -
                                          (problem.n_left + 1)]);
        j += 5;
      }
      for (int l = 0; l < left; l++) {
        obj_start -=
            (double)
                problem.crossing_matrix[solution[idx] - (problem.n_left + 1)]
                                       [solution[j + l] - (problem.n_left + 1)];
      }
    }
    solution.erase(solution.begin() + idx);
  }
  std::shuffle(new_nodes.begin(), new_nodes.end(), gen);
  return new_nodes;
}

std::vector<int> select_random(Problem &problem, std::vector<int> &solution,
                               int k, std::mt19937 &gen, double &obj_start) {
  std::uniform_int_distribution<int> dist;
  std::vector<int> new_nodes;
  int idx;
  while (new_nodes.size() < k) {
    dist = std::uniform_int_distribution<int>(0, solution.size() - 1);
    idx = dist(gen);

    new_nodes.push_back(solution[idx]);
    // substract crossings of the node to remove from objective value
    if (problem.graph[solution[idx] - problem.n_left - 1].size() > 0) {
      /*for(int i=0; i<idx; i++){
              //obj_start -=
      (double)problem.crossing_matrix[solution[i]-(problem.n_left+1)][solution[idx]-(problem.n_left+1)];
              obj_start -=
      (double)problem.crossing_matrix_transposed[solution[idx]-(problem.n_left+1)][solution[i]-(problem.n_left+1)];
      }
      for(int i=idx+1; i<solution.size(); i++){
              obj_start -=
      (double)problem.crossing_matrix[solution[idx]-(problem.n_left+1)][solution[i]-(problem.n_left+1)];
      }*/

      // unrolled version of the two previous loops
      int repeat = idx / 5;
      int left = idx % 5;
      int j = 0;
      while (repeat--) {
        obj_start -= (double)(problem.crossing_matrix_transposed
                                  [solution[idx] - (problem.n_left + 1)]
                                  [solution[j] - (problem.n_left + 1)] +
                              problem.crossing_matrix_transposed
                                  [solution[idx] - (problem.n_left + 1)]
                                  [solution[j + 1] - (problem.n_left + 1)] +
                              problem.crossing_matrix_transposed
                                  [solution[idx] - (problem.n_left + 1)]
                                  [solution[j + 2] - (problem.n_left + 1)] +
                              problem.crossing_matrix_transposed
                                  [solution[idx] - (problem.n_left + 1)]
                                  [solution[j + 3] - (problem.n_left + 1)] +
                              problem.crossing_matrix_transposed
                                  [solution[idx] - (problem.n_left + 1)]
                                  [solution[j + 4] - (problem.n_left + 1)]);
        j += 5;
      }
      for (int l = 0; l < left; l++) {
        obj_start -= (double)problem.crossing_matrix_transposed
                         [solution[idx] - (problem.n_left + 1)]
                         [solution[j + l] - (problem.n_left + 1)];
      }

      j = idx + 1;
      repeat = ((int)solution.size() - j) / 5;
      left = ((int)solution.size() - j) % 5;
      while (repeat--) {
        obj_start -=
            (double)(problem
                         .crossing_matrix[solution[idx] - (problem.n_left + 1)]
                                         [solution[j] - (problem.n_left + 1)] +
                     problem
                         .crossing_matrix[solution[idx] - (problem.n_left + 1)]
                                         [solution[j + 1] -
                                          (problem.n_left + 1)] +
                     problem
                         .crossing_matrix[solution[idx] - (problem.n_left + 1)]
                                         [solution[j + 2] -
                                          (problem.n_left + 1)] +
                     problem
                         .crossing_matrix[solution[idx] - (problem.n_left + 1)]
                                         [solution[j + 3] -
                                          (problem.n_left + 1)] +
                     problem
                         .crossing_matrix[solution[idx] - (problem.n_left + 1)]
                                         [solution[j + 4] -
                                          (problem.n_left + 1)]);
        j += 5;
      }
      for (int l = 0; l < left; l++) {
        obj_start -=
            (double)
                problem.crossing_matrix[solution[idx] - (problem.n_left + 1)]
                                       [solution[j + l] - (problem.n_left + 1)];
      }
    }
    solution.erase(solution.begin() + idx);
  }
  return new_nodes;
}

std::vector<int> lns(Problem &problem, int seed,
                     std::vector<int> start_solution, long n_iter,
                     int time_limit, bool verbose = true) {
  double elapsed;
  struct timespec end, start_time, current;
  clock_gettime(CLOCK_MONOTONIC, &start_time);

  std::mt19937 gen = std::mt19937(seed);
  std::uniform_int_distribution<int> dist_k_out =
      std::uniform_int_distribution<int>(100, 150);
  int k_out_block = std::min(50, (int)start_solution.size() - 1);
  std::uniform_real_distribution<double> dist_p =
      std::uniform_real_distribution<double>(0.0, 1.0);
  std::vector<int> solution_iter, new_nodes;
  std::vector<int> best_solution = start_solution;
  std::vector<int> best_solution_global = start_solution;
  double obj;
  double best_obj = problem.eval(best_solution);
  double best_obj_global = best_obj;
  std::uniform_int_distribution<int> dist;
  if (best_obj == 0.0) {
    return best_solution;
  }
  long iter = 0;
  int no_improvement = 0;
  while (true) {
    if (time_limit > 0) {
      clock_gettime(CLOCK_MONOTONIC, &current);
      elapsed = current.tv_sec - start_time.tv_sec +
                (current.tv_nsec - start_time.tv_nsec) * 1e-9;
      if (elapsed >= time_limit) {
        break;
      }
    }
    if (iter == n_iter) {
      break;
    }
    if (verbose && iter % 100 == 0) {
      printf("%ld %f (%f)\n", iter, best_obj, best_obj_global);
    }
    int k_out = std::min(dist_k_out(gen), (int)start_solution.size());
    solution_iter = best_solution;

    double p = dist_p(gen);
    double obj_start = best_obj;
    if (p < 0.2) {
      new_nodes = select_random(problem, solution_iter, k_out, gen, obj_start);
    } else {
      new_nodes =
          select_block(problem, solution_iter, k_out_block, gen, obj_start);
    }

    solution_iter =
        greedy_insertion(problem, solution_iter, new_nodes, obj, obj_start);

    if (obj < best_obj) {
      best_obj = obj;
      best_solution = solution_iter;
      if (best_obj < best_obj_global) {
        best_obj_global = best_obj;
        best_solution_global = solution_iter;
      }
      if (best_obj == 0.0) {
        break;
      }
      no_improvement = 0;
    } else {
      no_improvement += 1;
      if (no_improvement == 2000) {

        if (verbose == true) {
          printf("RESTART\n");
        }

        solution_iter = best_solution_global;
        obj_start = best_obj_global;

        new_nodes = select_random(problem, solution_iter,
                                  std::min(20, (int)solution_iter.size()), gen,
                                  obj_start);
        solution_iter = random_insertion(problem, solution_iter, new_nodes, obj,
                                         gen, obj_start);

        best_obj = obj;
        best_solution = solution_iter;

        no_improvement = 0;
      }
    }
    iter++;
  }
  if (verbose) {
    printf("Final: %f\n", best_obj);
  }
  return best_solution_global;
}

std::vector<int> solve(std::string problem_file, int seed, int time_limit,
                       bool verbose, bool read_problem_from_cin = false) {
  int threshold_segmentation =
      15000; // if number of nodes on right side exceeds this value, we compute
             // solution in segments
  int lns_iter_segmentation = 1000; // number of LNS iterations computed per
                                    // segment if segmentation is used
  double elapsed;
  struct timespec end, start_time, current;
  clock_gettime(CLOCK_MONOTONIC, &start_time);

  Problem problem;
  if (read_problem_from_cin) {
    problem.read_from_cin(threshold_segmentation);
  } else {
    problem.read_in(problem_file, threshold_segmentation);
  }
  double obj;
  std::vector<int> solution(problem.n_right);
  std::iota(solution.begin(), solution.end(), problem.n_left + 1);
  clock_gettime(CLOCK_MONOTONIC, &current);
  elapsed = current.tv_sec - start_time.tv_sec +
            (current.tv_nsec - start_time.tv_nsec) * 1e-9;
  if (elapsed >= time_limit) {

    return solution;
  }

  if (problem.n_right < threshold_segmentation) {
    // solution = barycenter(problem);

    std::vector<int> empty;
    solution = greedy_insertion(problem, empty, solution, obj);
    if (verbose) {
      printf("Greedy unsorted: %f\n", obj);
    }
    clock_gettime(CLOCK_MONOTONIC, &current);
    elapsed = current.tv_sec - start_time.tv_sec +
              (current.tv_nsec - start_time.tv_nsec) * 1e-9;
    if (elapsed >= time_limit) {
      return solution;
    }

    solution = lns(problem, seed, solution, std::numeric_limits<long>::max(),
                   time_limit - (int)elapsed, verbose);
  } else {
    // solution = barycenter(problem);
    if (problem.n_edges > 300000) {
      // if number of edges is high, the greedy insertion without crossing
      // matrix is too expensive
      solution = barycenter(problem);
    } else {
      std::vector<int> empty;
      solution = greedy_insertion(problem, empty, solution, obj);
      if (verbose) {
        printf("Greedy unsorted: %f\n", obj);
      }
    }

    clock_gettime(CLOCK_MONOTONIC, &current);
    elapsed = current.tv_sec - start_time.tv_sec +
              (current.tv_nsec - start_time.tv_nsec) * 1e-9;
    if (elapsed >= time_limit) {
      return solution;
    }
    int segment_size = 2000;
    int iter = 0;
    while (true) {
      if (verbose) {
        printf("Segment Size: %d\n", segment_size);
      }
      int num_segments = problem.n_right / segment_size + 1;
      for (int i = 0; i < num_segments; i++) {
        vector<int> subres;
        int start = std::max(
            0, std::min(i * segment_size, problem.n_right - segment_size - 1));
        int end = std::min((i + 1) * segment_size, problem.n_right);
        std::copy(solution.begin() + start, solution.begin() + end,
                  std::back_inserter(subres));

        Problem subproblem = problem.get_subproblem(subres);

        vector<int> subres_idx;
        for (int n = 0; n < subres.size(); n++) {
          subres_idx.push_back(problem.n_left + n + 1);
        }
        subproblem.compute_crossing_matrix();
        double obj_sub = subproblem.eval(subres_idx);
        if (verbose) {
          printf("Obj segment %d/%d: %f\n", i + 1, num_segments, obj_sub);
        }
        clock_gettime(CLOCK_MONOTONIC, &current);
        elapsed = current.tv_sec - start_time.tv_sec +
                  (current.tv_nsec - start_time.tv_nsec) * 1e-9;
        if (elapsed >= time_limit) {
          return solution;
        }
        if (obj_sub > 0) {
          vector<int> sol_lns =
              lns(subproblem, seed + iter, subres_idx, lns_iter_segmentation,
                  time_limit - (int)elapsed, verbose);
          for (int j = 0; j < sol_lns.size(); j++) {
            solution[start + j] = subres[sol_lns[j] - (problem.n_left + 1)];
          }
        }

        clock_gettime(CLOCK_MONOTONIC, &current);
        elapsed = current.tv_sec - start_time.tv_sec +
                  (current.tv_nsec - start_time.tv_nsec) * 1e-9;
        if (elapsed >= time_limit) {
          return solution;
        }
      }
      iter += 1;
      segment_size = std::min(segment_size + 1000, 5000);
    }
  }
  return solution;
}

int main(int argc, char **argv) {
  int PRODUCTION_STATE = true;
  int SEED = 1903;
  int TIME_LIMIT = 5 * 60 - 15;

  if (!PRODUCTION_STATE) {
    for (int problem_idx = 1; problem_idx <= 100; problem_idx++) {
      std::string problem_file =
          "../public_instances/" + std::to_string(problem_idx) + ".gr";
      Problem problem;
      problem.read_in(problem_file, 15000);

      printf("Problem: %d (%d, %d)\n", problem_idx, problem.n_left,
             problem.n_right);
      std::vector<int> solution;
      solution = solve(problem_file, SEED, TIME_LIMIT, true);
      printf("Obj: %f\n", problem.eval(solution));
    }
  } else {
    std::vector<int> solution;
    solution = solve("", SEED, TIME_LIMIT, false, true);
    for (int i = 0; i < solution.size(); i++) {
      printf("%d\n", solution[i]);
    }
  }
  return 0;
}
