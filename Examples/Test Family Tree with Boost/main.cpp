//=======================================================================
// Copyright 2001 Jeremy G. Siek, Andrew Lumsdaine, Lie-Quan Lee,
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <chrono>

using namespace std::chrono;
using namespace std;

enum family
{
    Glanz,
    Freya,
    Sierra,
    Zefir,
    Iyzumrae,
    Mithra,
    Solreya,
    Catenary,
    Znane,
    N
};
int main()
{
	using namespace boost;
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	const char* name[] = { "Glanz", "Freya", "Sierra", "Zefir", "Iyzumrae",
        "Mithra", "Solreya","Catenary", "Znane" };

	adjacency_list<> g(N);
	add_edge(Glanz, Zefir, g);
	add_edge(Glanz, Iyzumrae, g);
	add_edge(Glanz, Mithra, g);
	add_edge(Glanz, Solreya, g);
	add_edge(Glanz, Catenary, g);
	add_edge(Glanz, Znane, g);
	add_edge(Sierra, Znane, g);
	add_edge(Freya, Zefir, g);
	add_edge(Freya, Iyzumrae, g);
	add_edge(Freya, Mithra, g);
	add_edge(Freya, Solreya, g);
	add_edge(Freya, Catenary, g);
	
	graph_traits< adjacency_list<> >::vertex_iterator i, end;
	graph_traits< adjacency_list<> >::adjacency_iterator ai, a_end;
	auto index_map = get(vertex_index, g);

	for (boost::tie(i, end) = vertices(g); i != end; ++i)
	{
	cout << name[get(index_map, *i)];
	boost::tie(ai, a_end) = adjacent_vertices(*i, g);
	if (ai == a_end)
	{
		cout << " has no children";
	} 
	else
        {
		cout << " is the parent of ";
	}
	for (; ai != a_end; ++ai)
        {
		cout << name[get(index_map, *ai)];
		if (boost::next(ai) != a_end)
		{
			cout << ", ";
		}
	}
        cout << endl;
    	}

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return EXIT_SUCCESS;
}
