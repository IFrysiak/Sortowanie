#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <string>
#include <type_traits>
#include <cmath>
#include <fstream>

using namespace std;

// ARRAYS ----------------------------------------------------------------------------------------------------------------------
template <typename T>
T* generate_random_array(int size) {
    T* arr = new T[size];
    random_device rd;
    mt19937 gen(rd());

    if constexpr (is_integral_v<T>) {
        uniform_int_distribution<T> dist(1, size * 10);
        for (int i = 0; i < size; ++i) {
            arr[i] = dist(gen);
        }
    } else if constexpr (is_floating_point_v<T>) {
        uniform_real_distribution<T> dist(1.0, static_cast<T>(size * 10));
        for (int i = 0; i < size; ++i) {
            arr[i] = dist(gen);
        }
    } else if constexpr (is_same_v<T, string>) {
        const string chars = "abcdefghijklmnopqrstuvwxyz";
        uniform_int_distribution<int> char_dist(0, chars.size() - 1);
        uniform_int_distribution<int> len_dist(3, 10);
        for (int i = 0; i < size; ++i) {
            int str_len = len_dist(gen);
            string s;
            for (int j = 0; j < str_len; ++j) {
                s += chars[char_dist(gen)];
            }
            arr[i] = s;
        }
    }

    return arr;
}

template <typename T>
T* generate_partially_sorted_array(int size, double sorted_percentage) {
    T* arr = generate_random_array<T>(size);
    int sorted_elements = static_cast<int>(size * sorted_percentage / 100.0);
    sort(arr, arr + sorted_elements);
    return arr;
}

template <typename T>
T* generate_reverse_sorted_array(int size) {
    T* arr = generate_random_array<T>(size);
    sort(arr, arr + size, greater<T>());
    return arr;
}

template <typename T>
void print_array(const T* arr, int size, int limit = 20) {
    cout << "[";
    int to_print = min(size, limit);
    for (int i = 0; i < to_print; ++i) {
        cout << arr[i];
        if (i < to_print - 1) cout << ", ";
    }
    if (to_print < size) cout << ", ...";
    cout << "]" << endl;
}

template <typename T>
void delete_array(T* arr) {
    delete[] arr;
}

// SORTING ALGORTYHMS -----------------------------------------------------------------------------------------------------------------------------

// MERGE SORT

// Funkcja pomocnicza: scalanie dwóch podtablic
template <typename T>
void merge(T* S1, int n1, T* S2, int n2, T* S) {
    int i = 0, j = 0, k = 0;
    while (i < n1 && j < n2) {
        if (S1[i] <= S2[j])
            S[k++] = S1[i++];
        else
            S[k++] = S2[j++];
    }
    while (i < n1)
        S[k++] = S1[i++];
    while (j < n2)
        S[k++] = S2[j++];
}

// Funkcja główna: sortowanie przez scalanie
template <typename T>
void merge_sort(T* S, int n) {
    if (n <= 1) return; // już posortowana

    int n1 = n / 2;
    int n2 = n - n1;

    // Tworzymy dwie tymczasowe tablice
    T* S1 = new T[n1];
    T* S2 = new T[n2];

    // Kopiujemy dane
    for (int i = 0; i < n1; ++i)
        S1[i] = S[i];
    for (int i = 0; i < n2; ++i)
        S2[i] = S[n1 + i];

    // Rekurencyjnie sortujemy obie części
    merge_sort(S1, n1);
    merge_sort(S2, n2);

    // Scalanie dwóch posortowanych tablic do oryginalnej
    merge(S1, n1, S2, n2, S);

    // Usuwamy tymczasowe tablice
    delete[] S1;
    delete[] S2;
}
// QUICK SORT

// Funkcja pomocnicza: krok sortowania quicksort
template <typename T>
void quickSortStep(T* S, int a, int b) {
    if (a >= b) return; // 0 lub 1 element - nic nie trzeba robić

    T pivot = S[a + (b - a) / 2]; // pivot jako środkowy element

    int l = a;
    int r = b;

    while (l <= r) {
        while (S[l] < pivot) l++;
        while (S[r] > pivot) r--;
        if (l <= r) {
            swap(S[l], S[r]);
            l++;
            r--;
        }
    }

    // Rekurencyjnie sortujemy dwie części
    if (a < r) quickSortStep(S, a, r);
    if (l < b) quickSortStep(S, l, b);
}

// Funkcja główna: quicksort dla tablicy
template <typename T>
void quickSort(T* S, int n) {
    if (n <= 1) return; // już posortowana
    quickSortStep(S, 0, n - 1);
}

// HEAP SORT FOR INTROSORT
// HEAPSORT
template <typename T>
void heapify(T* arr, int n, int i) {
    int largest = i;
    int l = 2*i + 1;
    int r = 2*i + 2;

    if (l < n && arr[l] > arr[largest])
        largest = l;

    if (r < n && arr[r] > arr[largest])
        largest = r;

    if (largest != i) {
        swap(arr[i], arr[largest]);
        heapify(arr, n, largest);
    }
}

template <typename T>
void heapSort(T* arr, int n) {
    // Budujemy kopiec (heap)
    for (int i = n / 2 - 1; i >= 0; --i)
        heapify(arr, n, i);

    // Wyciągamy elementy z kopca jeden po drugim
    for (int i = n - 1; i > 0; --i) {
        swap(arr[0], arr[i]);
        heapify(arr, i, 0);
    }
}

// INTROSORT

// Funkcja pomocnicza do introsortu
template <typename T>
void introsortUtil(T* arr, int begin, int end, int depthLimit) {
    int size = end - begin + 1;
    if (size <= 16) {
        quickSortStep(arr, begin, end);
        return;
    }

    if (depthLimit == 0) {
        heapSort(arr + begin, size);
        return;
    }

    T pivot = arr[end];
    int l = begin;
    int r = end - 1;

    while (l <= r) {
        while (l <= r && !(pivot < arr[l])) l++;
        while (r >= l && !(arr[r] < pivot)) r--;
        if (l < r)
            swap(arr[l], arr[r]);
    }
    swap(arr[l], arr[end]);

    introsortUtil(arr, begin, l - 1, depthLimit - 1);
    introsortUtil(arr, l + 1, end, depthLimit - 1);
}

// Funkcja główna: introsort
template <typename T>
void introsort(T* arr, int n) {
    int depthLimit = 2 * log2(n);
    introsortUtil(arr, 0, n - 1, depthLimit);
}

// MESURE TIME FUNCTION
template<typename T, typename SortFunc>
double measure_sort_time(T* arr, int size, SortFunc sort_func) {
    auto start = chrono::high_resolution_clock::now();
    sort_func(arr, size);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> elapsed = end - start;
    return elapsed.count();
}

// MAIN -----------------------------------------------------------------------------------------------------------------------------
int main() {
    const int num_sizes = 9;
    int sizes[num_sizes] = {100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000};
    const int num_percentages = 6;
    double percentages[num_percentages] = {25, 50, 75, 95, 99, 99.7};
    const int repetitions = 100;

    ofstream outfile("results.csv");
    outfile << "Algorithm,ArrayType,ArraySize,OrderType,AverageTimeMs\n";

    for (int size_idx = 0; size_idx < num_sizes; ++size_idx) {
        int size = sizes[size_idx];
        
        // ==================== RANDOM =========================
        double merge_sort_total_time_random = 0.0;
        double quick_sort_total_time_random = 0.0;
        double intro_sort_total_time_random = 0.0;

        for (int rep = 0; rep < repetitions; ++rep) {
            int* random_arr = generate_random_array<int>(size);

            {
                int* arr_copy = new int[size];
                copy(random_arr, random_arr + size, arr_copy);
                double time_ms = measure_sort_time(arr_copy, size, merge_sort<int>);
                merge_sort_total_time_random += time_ms;
                delete[] arr_copy;
            }
            {
                int* arr_copy = new int[size];
                copy(random_arr, random_arr + size, arr_copy);
                double time_ms = measure_sort_time(arr_copy, size, quickSort<int>);
                quick_sort_total_time_random += time_ms;
                delete[] arr_copy;
            }
            {
                int* arr_copy = new int[size];
                copy(random_arr, random_arr + size, arr_copy);
                double time_ms = measure_sort_time(arr_copy, size, introsort<int>);
                intro_sort_total_time_random += time_ms;
                delete[] arr_copy;
            }

            delete_array(random_arr);
        }

        outfile << "MergeSort,Random," << size << ",Random," << (merge_sort_total_time_random / repetitions) << "\n";
        outfile << "QuickSort,Random," << size << ",Random," << (quick_sort_total_time_random / repetitions) << "\n";
        outfile << "IntroSort,Random," << size << ",Random," << (intro_sort_total_time_random / repetitions) << "\n";
        
        /*
        // ==================== PARTIALLY SORTED =========================
        for (int pct_idx = 0; pct_idx < num_percentages; ++pct_idx) {
            double percent = percentages[pct_idx];

            double merge_sort_total_time_partial = 0.0;
            double quick_sort_total_time_partial = 0.0;
            double intro_sort_total_time_partial = 0.0;

            for (int rep = 0; rep < repetitions; ++rep) {
                int* partial_arr = generate_partially_sorted_array<int>(size, percent);

                {
                    int* arr_copy = new int[size];
                    copy(partial_arr, partial_arr + size, arr_copy);
                    double time_ms = measure_sort_time(arr_copy, size, merge_sort<int>);
                    merge_sort_total_time_partial += time_ms;
                    delete[] arr_copy;
                }
                {
                    int* arr_copy = new int[size];
                    copy(partial_arr, partial_arr + size, arr_copy);
                    double time_ms = measure_sort_time(arr_copy, size, quickSort<int>);
                    quick_sort_total_time_partial += time_ms;
                    delete[] arr_copy;
                }
                {
                    int* arr_copy = new int[size];
                    copy(partial_arr, partial_arr + size, arr_copy);
                    double time_ms = measure_sort_time(arr_copy, size, introsort<int>);
                    intro_sort_total_time_partial += time_ms;
                    delete[] arr_copy;
                }

                delete_array(partial_arr);
            }

            outfile << "MergeSort,PartialSorted," << size << "," << percent << "%," << (merge_sort_total_time_partial / repetitions) << "\n";
            outfile << "QuickSort,PartialSorted," << size << "," << percent << "%," << (quick_sort_total_time_partial / repetitions) << "\n";
            outfile << "IntroSort,PartialSorted," << size << "," << percent << "%," << (intro_sort_total_time_partial / repetitions) << "\n";
        }
        */
        /*
        // ==================== REVERSE SORTED =========================
        double merge_sort_total_time_reverse = 0.0;
        double quick_sort_total_time_reverse = 0.0;
        double intro_sort_total_time_reverse = 0.0;

        for (int rep = 0; rep < repetitions; ++rep) {
            int* reverse_arr = generate_reverse_sorted_array<int>(size);

            {
                int* arr_copy = new int[size];
                copy(reverse_arr, reverse_arr + size, arr_copy);
                double time_ms = measure_sort_time(arr_copy, size, merge_sort<int>);
                merge_sort_total_time_reverse += time_ms;
                delete[] arr_copy;
            }
            {
                int* arr_copy = new int[size];
                copy(reverse_arr, reverse_arr + size, arr_copy);
                double time_ms = measure_sort_time(arr_copy, size, quickSort<int>);
                quick_sort_total_time_reverse += time_ms;
                delete[] arr_copy;
            }
            {
                int* arr_copy = new int[size];
                copy(reverse_arr, reverse_arr + size, arr_copy);
                double time_ms = measure_sort_time(arr_copy, size, introsort<int>);
                intro_sort_total_time_reverse += time_ms;
                delete[] arr_copy;
            }

            delete_array(reverse_arr);
        }

        outfile << "MergeSort,ReverseSorted," << size << ",Reverse," << (merge_sort_total_time_reverse / repetitions) << "\n";
        outfile << "QuickSort,ReverseSorted," << size << ",Reverse," << (quick_sort_total_time_reverse / repetitions) << "\n";
        outfile << "IntroSort,ReverseSorted," << size << ",Reverse," << (intro_sort_total_time_reverse / repetitions) << "\n";
        */   
    }

    outfile.close();
    cout << "Eksperymenty zakończone. Wyniki zapisane w pliku 'results.csv'.\n";

    return 0;
}
// Kompilacja: g++ -std=c++2a test.cpp
// Uruchomienie: .\a.exe