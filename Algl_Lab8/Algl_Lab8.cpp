#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <random>
#include <numeric>

using namespace std;

/**
 * A class that represents a MaxHeap data structure.
 */
class MaxHeap {
private:
    std::vector<unsigned int> heap;

    // Helper functions

    /**
     * Helper function to calculate the index of a node's parent node.
     * @param i The index of the node for which to calculate the parent index.
     * @return The index of the parent node.
     */
    int parent(int i) { return (i - 1) / 2; }

    /**
     * Helper function to calculate the index of a node's left child node.
     * @param i The index of the node for which to calculate the left child index.
     * @return The index of the left child node.
     */
    int leftChild(int i) { return 2 * i + 1; }

    /**
     * Helper function to calculate the index of a node's right child node.
     * @param i The index of the node for which to calculate the right child index.
     * @return The index of the right child node.
     */
    int rightChild(int i) { return 2 * i + 2; }

    /**
     * Moves an element up the heap to its correct position.
     * @param i The index of the element to be moved up the heap.
     */
    void heapifyUp(int i) {
        while (i != 0 && heap[i] > heap[parent(i)]) {
            swap(heap[i], heap[parent(i)]);
            i = parent(i);
        }
    }

    /**
     * Moves an element down the heap to its correct position.
     * @param i The index of the element to be moved down the heap.
     */
    void heapifyDown(int i) {
        int maxIndex = i;

        int l = leftChild(i);
        if (l < heap.size() && heap[l] > heap[maxIndex]) {
            maxIndex = l;
        }

        int r = rightChild(i);
        if (r < heap.size() && heap[r] > heap[maxIndex]) {
            maxIndex = r;
        }

        if (i != maxIndex) {
            swap(heap[i], heap[maxIndex]);
            heapifyDown(maxIndex);
        }
    }

public:
    /**
    * Constructor that creates an empty MaxHeap.
    */
    MaxHeap() {}

    /**
     * Inserts an element into the MaxHeap.
     * @param val The value of the element to be inserted.
     */
    void insert(int val) {
        heap.push_back(val);
        heapifyUp(heap.size() - 1);
    }

    /**
     * Removes the maximum element from the MaxHeap.
     * @return The value of the maximum element that was removed.
     */
    int deleteMax() {
        if (heap.empty()) {
            throw runtime_error("Heap is empty");
        }
        int max = heap[0];
        heap[0] = heap.back();
        heap.pop_back();
        if (!heap.empty()) {
            heapifyDown(0);
        }
        return max;
    }

    /**
     * Returns the maximum element in the MaxHeap.
     * @return The value of the maximum element.
     */
    int findMax() {
        if (heap.size() == 0) {
            throw runtime_error("Heap is empty");
        }
        return heap[0];
    }

    /**
     * Returns the number of elements in the MaxHeap.
     * @return The number of elements in the MaxHeap.
     */
    int size() const {
        return heap.size();
    }
};


// Node structure for binomial heap
struct Node {
    int val; // value of node
    int degree; // degree of node
    Node* parent; // parent of node
    Node* child; // child of node
    Node* sibling; // sibling of node

    Node(int val, int degree, Node* parent, Node* child, Node* sibling) : val(val), degree(degree), parent(parent), child(child), sibling(sibling) {}
};

class BinomialMaxHeap {
private:
    struct Node {
        int val;
        int degree;
        Node* parent;
        Node* child;
        Node* sibling;

        Node(int val) : val(val), degree(0), parent(nullptr), child(nullptr), sibling(nullptr) {}
    };

    vector<Node*> trees;
    int size;

    /**
     * Merge two binomial trees of the same degree.
     *
     * @param t1 The first binomial tree to merge.
     * @param t2 The second binomial tree to merge.
     * @return The root node of the merged binomial tree.
     */
    Node* mergeTrees(Node* t1, Node* t2) {
        if (t1->val < t2->val) swap(t1, t2);
        t2->parent = t1;
        t2->sibling = t1->child;
        t1->child = t2;
        t1->degree++;
        return t1;
    }

    /**
     * Merge two binomial heaps.
     *
     * @param h1 The first binomial heap to merge.
     * @param h2 The second binomial heap to merge.
     * @return The merged binomial heap.
     */
    BinomialMaxHeap mergeHeaps(BinomialMaxHeap h1, BinomialMaxHeap h2) {
        BinomialMaxHeap h;
        h.trees.reserve(max(h1.trees.size(), h2.trees.size()) + 1);

        Node* carry = nullptr;
        auto mergeTreesToHeap = [&](Node* t) {
            if (carry != nullptr) {
                h.trees.push_back(carry);
                carry = nullptr;
            }
            while (!h.trees.empty() && h.trees.back()->degree == t->degree) {
                t = mergeTrees(h.trees.back(), t);
                h.trees.pop_back();
            }
            h.trees.push_back(t);
        };

        auto mergeSiblingTreesToHeap = [&](Node* t) {
            while (t != nullptr) {
                Node* sib = t->sibling;
                t->parent = nullptr;
                t->sibling = nullptr;
                mergeTreesToHeap(t);
                t = sib;
            }
        };

        if (!h1.trees.empty() && !h2.trees.empty()) {
            mergeSiblingTreesToHeap(h1.trees.front());
            mergeSiblingTreesToHeap(h2.trees.front());
        }
        else if (!h1.trees.empty()) {
            mergeSiblingTreesToHeap(h1.trees.front());
        }
        else if (!h2.trees.empty()) {
            mergeSiblingTreesToHeap(h2.trees.front());
        }

        if (carry != nullptr) {
            h.trees.push_back(carry);
        }

        return h;
    }

public:
    /**
     * Constructor for creating an empty binomial heap.
     */
    BinomialMaxHeap() : size(0) {}

    /**
     * Inserts a value into the binomial heap.
     *
     * @param val The value to insert.
     */
    void insert(int val) {
        BinomialMaxHeap h;
        h.trees.push_back(new Node(val));
        *this = mergeHeaps(*this, h);
        size++;
    }

    /**
     * Finds the maximum value in the binomial heap.
     *
     * @return The maximum value in the binomial heap.
     */
    int findMax() const {
        return (*max_element(trees.begin(), trees.end(), [](Node* a, Node* b) { return a->val < b->val; }))->val;
    }
    
    /**
    * Deletes the maximum value from the binomial heap.
    */
    void deleteMax() {
        if (size == 0) return;

        // Find the tree with the maximum value
        auto maxTreeIt = max_element(trees.begin(), trees.end(), [](Node* a, Node* b) { return a->val > b->val; });
        Node* maxTree = *maxTreeIt;

        // Remove the max tree from the heap
        trees.erase(maxTreeIt);
        size -= (1 << maxTree->degree);

        // Create a new heap from the max tree's children
        BinomialMaxHeap h;
        if (maxTree->child != nullptr) {
            h.trees.resize(maxTree->degree);
            int i = 0;
            for (Node* child = maxTree->child; child != nullptr; child = child->sibling) {
                child->parent = nullptr;
                h.trees[i++] = child;
            }
            *this = mergeHeaps(*this, h);
        }

        // Delete the max tree
        delete maxTree;
    }

};


std::vector<int> genHeapVector(int N) {
    std::vector<int> vector(pow(10, N));
    mt19937 engine(time(nullptr));
    uniform_int_distribution<int> gen(0, pow(10, N));
    for (auto& el : vector) {
        el = gen(engine);
    }

    return vector;
}

std::vector<int> genVector(int N) {
    std::vector<int> vector(1000);
    mt19937 engine(time(nullptr));
    uniform_int_distribution<int> gen(0, pow(10, N));
    for (auto& el : vector) {
        el = gen(engine);
    }

    return vector;
}


int main() {
    ofstream fout("delete.txt");

    MaxHeap maxh;
    std::vector<int> vectorrr(1000);
    std::iota(vectorrr.begin(), vectorrr.end(), 1);
    for (auto& el : vectorrr)
    {
        maxh.insert(el);
    }
    std::cout << maxh.findMax() << endl;
    std::cout << maxh.deleteMax() << endl;
    std::cout << maxh.findMax() << endl;
    std::cout << maxh.deleteMax() << endl;
    
    std::cout << "-------------------------" << endl;
    BinomialMaxHeap maxb;
    for (auto& el : vectorrr)
    {
        maxb.insert(el);
    }
    std::cout << maxb.findMax() << endl;
    maxb.deleteMax();
    std::cout << maxb.findMax() << endl;
    maxb.deleteMax();
    std::cout << maxb.findMax() << endl;

    for (int i = 3; i < 8; ++i) {
        BinomialMaxHeap binomialHeap;
        MaxHeap maxHeap;

        std::vector<int> vec = genHeapVector(i);
        for (auto el : vec) {
            binomialHeap.insert(el);
            maxHeap.insert(el);
        }

        /*time of inserting 1000 elements into maxHeap*/
        std::vector<double> findMaxTime6;
        std::vector<int> vec2 = genVector(i);
        for (auto& el1 : vec2) {
            chrono::high_resolution_clock::time_point start1 = chrono::high_resolution_clock::now();
            maxHeap.insert(el1);
            chrono::high_resolution_clock::time_point end1 = chrono::high_resolution_clock::now();
            chrono::duration<double> sec_diff1 = end1 - start1;
            findMaxTime6.push_back(sec_diff1.count());
        }
        double time6 = std::accumulate(findMaxTime6.begin(), findMaxTime6.end(), 0.0);
        double max6 = *std::max_element(findMaxTime6.begin(), findMaxTime6.end());
        //fout << "time inserting 1000 elements in maxHeap with " << pow(10, i) << " elements: " << time6 << endl;
        //fout << "worst time inserting 1 element in maxHeap with " << pow(10, i) << " elements: " << max6 << endl;

        /*time of inserting 1000 elements into binomialHeap*/
        std::vector<double> findMaxTime5;
        for (auto& el : vec2) {
            chrono::high_resolution_clock::time_point start2 = chrono::high_resolution_clock::now();
            binomialHeap.insert(el);
            chrono::high_resolution_clock::time_point end2 = chrono::high_resolution_clock::now();
            chrono::duration<double> sec_diff2 = end2 - start2;
            findMaxTime5.push_back(sec_diff2.count());
        }
        double time5 = std::accumulate(findMaxTime5.begin(), findMaxTime5.end(), 0.0);
        double max5 = *std::max_element(findMaxTime5.begin(), findMaxTime5.end());
        //fout << "time inserting 1000 elements in binomialHeap with " << pow(10, i) << " elements: " << time5 << endl;
        //fout << "worst time inserting 1 element in binomialHeap with " << pow(10, i) << " elements: " << max5 << endl;

        ///*time of 1000 findMax() operations in maxHeap*/
        //std::vector<double> findMaxTime;
        //for (int k = 0; k < 1000; ++k) {
        //    chrono::high_resolution_clock::time_point start3 = chrono::high_resolution_clock::now();
        //    maxHeap.findMax();
        //    chrono::high_resolution_clock::time_point end3 = chrono::high_resolution_clock::now();
        //    chrono::duration<double> sec_diff3 = end3 - start3;
        //    findMaxTime.push_back(sec_diff3.count());
        //}
        //double time1 = std::accumulate(findMaxTime.begin(), findMaxTime.end(), 0.0);
        //double max1 = *std::max_element(findMaxTime.begin(), findMaxTime.end());
        ///*fout << "time of 1000 findMax() operations in maxHeap with " << pow(10, i) << " elements: " << time1 << endl;
        //fout << "time of 1 worst findMax() operation in maxHeap with " << pow(10, i) << " elements: " << max1 << endl;*/


        ///*time of 1000 findMax() operations in binomialHeap*/
        //std::vector<double> findMaxTime2;
        //for (int k = 0; k < 1000; ++k) {
        //    chrono::high_resolution_clock::time_point start3 = chrono::high_resolution_clock::now();
        //    binomialHeap.findMax();
        //    chrono::high_resolution_clock::time_point end3 = chrono::high_resolution_clock::now();
        //    chrono::duration<double> sec_diff3 = end3 - start3;
        //    findMaxTime2.push_back(sec_diff3.count());
        //}
        //double time2 = std::accumulate(findMaxTime2.begin(), findMaxTime2.end(), 0.0);
        //double max2 = *std::max_element(findMaxTime2.begin(), findMaxTime2.end());
        ///*fout << "time of 1000 findMax() operations in in binomialHeap with " << pow(10, i) << " elements: " << time2 << endl;
        //fout << "time of 1 worst findMax() operation in binomialHeap with " << pow(10, i) << " elements: " << max2 << endl;*/

        /*time of 1000 deleteMax() operations in maxHeap*/
        std::vector<double> findMaxTime3;
        for (int k = 0; k < 1000; ++k) {
            chrono::high_resolution_clock::time_point start3 = chrono::high_resolution_clock::now();
            maxHeap.deleteMax();
            chrono::high_resolution_clock::time_point end3 = chrono::high_resolution_clock::now();
            chrono::duration<double> sec_diff3 = end3 - start3;
            findMaxTime3.push_back(sec_diff3.count());
        }
        double time3 = std::accumulate(findMaxTime3.begin(), findMaxTime3.end(), 0.0);
        double max3 = *std::max_element(findMaxTime3.begin(), findMaxTime3.end());
        fout << "time of 1000 deleteMax() operations in in maxHeap with " << pow(10, i) << " elements: " << time3 << endl;
        fout << "time of 1 worst deleteMax() operation in maxHeap with " << pow(10, i) << " elements: " << max3 << endl;

        /*time of 1000 deleteMax() operations in binomialHeap*/
        std::vector<double> findMaxTime4;
        for (int k = 0; k < 1000; ++k) {
            chrono::high_resolution_clock::time_point start3 = chrono::high_resolution_clock::now();
            binomialHeap.deleteMax();
            chrono::high_resolution_clock::time_point end3 = chrono::high_resolution_clock::now();
            chrono::duration<double> sec_diff3 = end3 - start3;
            findMaxTime4.push_back(sec_diff3.count());
        }
        double time4 = std::accumulate(findMaxTime4.begin(), findMaxTime4.end(), 0.0);
        double max4 = *std::max_element(findMaxTime4.begin(), findMaxTime4.end());
        fout << "time of 1000 deleteMax() operations in in binomialHeap with " << pow(10, i) << " elements: " << time4 << endl;
        fout << "time of 1 worst deleteMax() operation in binomialHeap with " << pow(10, i) << " elements: " << max4 << endl;


    }

    
    fout.close();
    return EXIT_SUCCESS;
}