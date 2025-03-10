#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <random>
#include <stack>
#include <filesystem>
#include <chrono>

using namespace std;
using namespace std::filesystem;

const long double PI = 3.14159265359, INF = 1e8;
int W;

double global_ratio;

struct Genom
{
    vector <bool> B;
    int weight;
    int value;

    Genom () : weight(0), value(0) {};
    Genom (vector <pair <int, int> > &a, unsigned int n, unsigned int i)
    {
        B.resize(n, 0);
        B[i] = 1;
        weight = a[i].second;
        value = a[i].first;
    }

    bool operator [](unsigned int i)
    {
        return B[i];
    }

    unsigned int size()
    {
        return B.size();
    }

    double get_ratio()
    {
        double ans = value;
        ans /= weight;
        return ans;
    }

    Genom* sum(vector <pair <int, int> > &a, Genom* g)
    {
        Genom* new_Genom = new Genom();
        new_Genom->B.resize(g->size());
        for (int i = 0; i < g->size(); i++)
        {
            int flag = (B[i] > 0 || g->B[i] > 0);
            new_Genom->B[i] = flag;
            new_Genom->weight += a[i].second * flag;
            new_Genom->value += a[i].first * flag;
        }
        return new_Genom;
    }
};

bool cmp(Genom* a, Genom* b)
{
//    return a->get_ratio() > b->get_ratio();
    return a->value > b->value;
}

struct Population
{
    vector <Genom*> G;

    Population (){};
    Population (vector <pair <int, int> > &a, unsigned int n, unsigned int max_P)
    {
        for (int i = 0; i < max_P; i++)
        {
            Genom* new_genom = new Genom(a, n, i%n);
            G.push_back(new_genom);
        }
    }
    ~Population ()
    {
        for (int i = 0; i < G.size(); i++)
        {
            delete G[i];
        }
    }

    Genom* operator [](unsigned int i)
    {
        return G[i];
    }

    void add(Genom* a)
    {
        G.push_back(a);
    }

    unsigned int size()
    {
        return G.size();
    }

    void gen_sort()
    {
        sort(G.begin(), G.end(), cmp);
    }

    void selection(unsigned int max_P)
    {
        vector <Genom*> new_G(max_P);
        for (int i = 0; i < max_P; i++)
        {
            new_G[i] = G[i];
        }
        G = new_G;
    }
};

int epoch(vector <pair <int, int> > &a, Population &P, unsigned int max_P)
{
//    cout << "start" << endl;
    for (int i = 0; i < max_P; i++)
    {
        for (int j = i+1; j < max_P; j++)
        {
            Genom* new_Genom = P[i]->sum(a, P[j]);
            if (new_Genom->weight < W)
                P.add(new_Genom);
        }
    }
//    for (int i = 0; i < P.size(); i++)
//    {
//        for (int j = 0; j < P[i]->size(); j++)
//        {
//            cout << (*P[i])[j] << " ";
//        }
//        cout << endl;
//    }

    P.gen_sort();

    P.selection(max_P);

    int value = 0;
    for (int i = 0; i < max_P; i++)
    {
        value = max(P[i]->value, value);
    }
    return value;
}

int main()
{
    W = 0;
    ifstream in("data/ks_19_0");

    int n, w;
    vector <pair <int, int> > a;

    double max_V = 0, max_W = 0;

    if (in.is_open())
    {
        in >> n >> w;
        W = w;
//        cout << n << " " << w << endl;

        for (int i = 0; i < n; i++)
        {
            int vi, wi;
            in >> vi >> wi;
            max_V += vi;
            max_W += wi;
            a.push_back(make_pair(vi, wi));
//            cout << vi << " " << wi << endl;
        }

//        cout << endl;
    }
    in.close();

///--------------------------

    global_ratio = max_V / max_W;

    int prev_value = 0;
    int new_value = -1;
    int flag = 0;

    unsigned int max_P = 4*n;
    Population P(a, n, max_P);

//    int value = 0;
//    for (int i = 0; i < max_P; i++)
//    {
//        value = max(P[i]->value, value);
//    }
//    cout << value << endl;

    while(flag < 2)
    {
        new_value = epoch(a, P, max_P);
        cout << new_value << endl;
//        for (int i = 0; i < max_P; i++)
//        {
//            for (int j = 0; j < P[i]->size(); j++)
//            {
//                cout << (*P[i])[j] << " ";
//            }
//            cout << endl;
//        }

        if (prev_value - new_value == 0) flag++;
        else flag = 0;

        prev_value = new_value;
    }

    cout << new_value << endl;

    return 0;
}
