/*
 * !/usr/bin/python
 * -*- coding: utf-8 -*-
 *
 * This file is part of pyunicorn.
 * Copyright (C) 2008--2019 Jonathan F. Donges and pyunicorn authors
 * URL: <http://www.pik-potsdam.de/members/donges/software>
 * License: BSD (3-clause)
 */

#include <math.h>
#define ARR_SIZE(arr) ( sizeof((arr)) / sizeof((arr[0])) )  // Always useful macro to get the array size!
// geo_network ================================================================


float *distance(float *D, int node, int *list_of_neighbors, int numNeighbors, int N ){
	float *d_node_all = malloc(sizeof( *d_node_all) * numNeighbors) ;
	for (int i=0; i<numNeighbors; i++) {
		d_node_all[i]=D[node*N + list_of_neighbors[i]];
		//		printf("d_node_all %d - %d: %2f \n", node, list_of_neighbors[i], d_node_all[i]);
	}
	return d_node_all;
}

// A utility function to swap to integers
void swap (int *a, int *b)
{
	int temp = *a;
	*a = *b;
	*b = temp;
}

// A utility function to print an array
void printArray (int arr[], int n)
{
	printf("Array Print: \n");
	for (int i = 0; i < n; i++){
		printf(" %d \n", arr[i]);
	}
	printf("\n");
}
// A utility function to print an array
void printFloatArray (float arr[], int n)
{
	printf("Array Print: \n");
	for (int i = 0; i < n; i++){
		printf(" %2f \n", arr[i]);
	}
	printf("\n");
}



// A function to generate a random permutation of neighbor_lists
void fisher_yates_shuffeling(int *arr, int n){
	srand (time(NULL));
	for (int i = n-1; i > 0; i--){
		// Start from the last element and swap
		// one by one. We don't need to run for
		// the first element that's why i > 0
		for (int i = n - 1; i > 0; i--)
		{
			// Pick a random index from 0 to i
			int j = rand() % (i + 1);
			// Swap arr[i] with the element at random index
			swap(&arr[i], &arr[j]);
		}
	}
}

int in_array(const int store[], const int storeSize, const int query) {
	for (int i=0; i<storeSize; ++i) {
		if (store[i] == query) {
			return i;
		}
	}
	return -1;
}

void _geo_model_1_fast(int iterations, float tolerance,
		short *A, float *D, int *link_list, int N, int E)
{

	//  Initialize random number generator
	srand48(time(0));
	int i, j, k=0, l=-1;
	int *list_all_neighbors=malloc(sizeof(*list_all_neighbors)*N) ;
	for (int tmp=0; tmp<N; tmp++) {
		list_all_neighbors[tmp]=tmp;
	}

	const int dim_link_list=2;
	//	for(int link=0; link<2*dim_list*E; link+=2*dim_list){
	//		printf("Link: %d, nodes: %d %d \n", link/(2*dim_list), link_list[link], link_list[link+dim_list]);
	//	}
	//	for(int i=0;i<N;i++){
	//		for(int j=0; j<N;j++){
	//			printf("%d ", A[i*N+j]);
	//		}
	//		printf("\n");
	//	}
	//	for(int i=0;i<N;i++){
	//		for(int j=0; j<N;j++){
	//			printf("%1f ", D[i*N+j]);
	//		}
	//		printf("\n");
	//	}
	//	return ;

	int q = 0;
	for (int u=0; u< iterations; u++) {

		int cond = 0;
		while (cond<1) {
			printf("cond: %d, q: %d \n", cond, q);
			q += 1;

			int rand_edge = floor(drand48() * E) ;
			int first_link_index=rand_edge*dim_link_list*2;
			i= link_list[first_link_index];
			j= link_list[first_link_index + dim_link_list];

			//			printf("Length link list: %ld, %d \n", ARR_SIZE(link_list), E);
			//			printf("Link index: %d,  Link i: %d, Link j: %d \n", first_link_index, i ,j );

			float *d_i_all;
			d_i_all=distance(D, i,list_all_neighbors, N, N );
			//			for (int it=0; it<N; it++) {
			//				printf("d_i_all %d - %d: %2f \n ", i, list_all_neighbors[it], d_i_all[it]);
			//			}

			printf("d_i_all\n");
			printFloatArray(d_i_all, N);

			// Create mask and provide distances in mask
			int nb_count=0;
			int mask[N];
			for (int d=0; d<N; d++){
				float Dist_j = d_i_all[d] - d_i_all[j];
				if (fabs(Dist_j) < tolerance * d_i_all[j]){
					mask[d]=1;
					nb_count++;
					//					printf("Node %d, nb_count %d, Dist_j: %2f \n", d, nb_count, Dist_j);
				} else {
					mask[d]=0;
				}
			}

			free(d_i_all);

			// We have to exclude the link i<->j from mask! Corresponds to False
			if (nb_count>2){
				mask[i] = 0;
				mask[j] = 0;
			} else {
				continue;
			}

			// Create list of possible neighbors
			int possible_nbs [nb_count-2];
			int tmp=0;
			for (int d=0; d<N; d++){
				if (mask[d]==1){
					possible_nbs[tmp]=d;
					tmp++;
				}
			}

			// Permutation of possible_nbs list
			if (nb_count>2){
				fisher_yates_shuffeling(possible_nbs,nb_count-2);
				printArray(possible_nbs, nb_count-2);
			}

			// Find neighbor in list of possible nbs with correct distance for nb_count -2 times as i, j are not counted)
			l = -1;
			for (int rk=0; rk<nb_count -2 ; rk++) {
				int nk_count=0;
				k=possible_nbs[rk]; // check now the neighbors of k
				int nbs_of_k[N]; // N is the maximum number of neighbor nodes in fully connected network
				for (int link=0; link<2*dim_link_list*E; link+=2*dim_link_list) {
					if (link_list[link]==k ) {
						nbs_of_k[nk_count]=link_list[link+dim_link_list];
						nk_count++;
					}
					// This if is necessary as link list is ordered from smaller node number to larger node number!
					if (link_list[link+dim_link_list]==k ) {
						nbs_of_k[nk_count]=link_list[link];
						nk_count++;
					}
				}

				printf("nbs_of_k k=%d \n", k);
				printArray(nbs_of_k, nk_count);

				// Test if possible neighbor occurs in nbs_of_k as well
				if ((in_array(nbs_of_k, nk_count, i)>-1) | (nk_count == 0) ) {
					//					printf("In Array test link i: %d or nk_count: %d, rk %d , while loop %d, step %d !\n", i, nk_count, rk, q, u);
					continue;
				} else {
					float *d_k_all;
					d_k_all= distance(D,k, nbs_of_k, nk_count,N);
					printf("d_k_all: k %d nk_count %d \n", k, nk_count);

					float *d_j_all;
					d_j_all= distance(D,j, nbs_of_k, nk_count,N);
					printf("d_j_all: j	 %d nk_count %d \n ", j, nk_count);

					int mask2[nk_count];
					int any_candidate=0;
					for(int d=0; d<nk_count;d++) {
						float Dist_k_j=d_k_all[d] - d_j_all[d];
						if (fabs(Dist_k_j)<tolerance*d_k_all[d]) {
							mask2[d]=1;
							any_candidate+=1;
						} else {
							mask2[d]=0;
						}
					}

					// Analyze list of possible candidates further
					if (any_candidate > 0) {
						//						printf("Any candidate: %d \n", any_candidate);
						int possible_candidates[any_candidate];
						int candidate_count=0;
						for(int tmp=0; tmp<nk_count; tmp++) {
							if (mask2[tmp]==1) {
								possible_candidates[candidate_count]=nbs_of_k[tmp];
								candidate_count++;
							}
						}

						printf("Possible_candidates: %d \n", candidate_count);
						printArray(possible_candidates, candidate_count);

						// Now check for neighbors of l_candidate
						int l_candidate=possible_candidates[rand() % any_candidate];
						int nl_count=0;
						int nbs_of_l[N];
						for(int link=0; link<2*dim_link_list*E; link+=2*dim_link_list){
							if (link_list[link]==l_candidate)
							{
								nbs_of_l[nl_count]=link_list[link+dim_link_list];
								nl_count++;
							}
							if (link_list[link+dim_link_list]==l_candidate)
							{
								nbs_of_l[nl_count]=link_list[link];
								nl_count++;
							}
						}
						printf("Possible_nl_list: %d \n", nl_count);
						printArray(nbs_of_l, nl_count);

						for (int nl=0; nl<nl_count; nl++) {
							// check if j is not in nbs_of_l
							if (in_array(nbs_of_l, nl_count, j)==-1) {
								l=l_candidate;
								printf("l swap found: l %d j %d \n", l,j);
								rk=nb_count; //  An exchange neighbor is found, leave loop over list of possible neighbors
								//								exit(0);
								break;
							}
						}
					}
					free(d_k_all);
					free(d_j_all);
				}
			}
			if (l==-1){
				printf("Return to beginning of while loop! \n");
				continue;  // Returns to beginning of while loop
			} else // Apply change
			{
				printf("Apply changes in adjacency matrix: i %d j %d; k %d l%d \n", i, j, k, l);
				// start changing links in adjacency matrix and link list!
				cond = 1;

				A[i*N + j] =  A[j*N + i] = 0;  // Delete link i<->j
				A[k*N + k] =  A[l*N + k] = 0;  // Delete link k<->l
				A[i*N + k] =  A[k*N + i] = 1;  // Add link i<->k
				A[j*N + l] =  A[l*N + j] = 1;  // Add link j<->l

				printf("Changed Adj matrix. \n");

				// Now find id of second_link_index k<->l
				int second_link_index=0;

				//	for(int link=0; link<2*dim_list*E; link+=2*dim_list){
				//		printf("Link: %d, nodes: %d %d \n", link/(2*dim_list), link_list[link], link_list[link+dim_list]);
				//	}

				for (int link=0; link<2*dim_link_list*E; link+=2*dim_link_list){
					if ( (link_list[link]==k && link_list[link+dim_link_list]==l ) ||  (link_list[link]==l && link_list[link+dim_link_list]==k ) ) {
						second_link_index=link;
						break;
					}
				}
				printf("Second link index found: %d \n", second_link_index);

				// Now update the link list to i<->k j<->l
				link_list[first_link_index] = i;
				link_list[first_link_index+ dim_link_list] = k;
				link_list[second_link_index] = j;
				link_list[second_link_index+dim_link_list] = l;

				printf("link list updated! \n") ;
			}
		}
	}
	printf("Trials %d, Iterations %d \n", q, iterations);
}


void _geo_model_2_fast(int iterations, float tolerance,
		short *A, float *D, int *link_list, int N, int E)
{

	//  Initialize random number generator
	srand48(time(0));
	int i, j, k=0, l=-1;
	int *list_all_neighbors=malloc(sizeof(*list_all_neighbors)*N) ;
	for (int tmp=0; tmp<N; tmp++) {
		list_all_neighbors[tmp]=tmp;
	}

	const int dim_link_list=2;
	//	for(int link=0; link<2*dim_list*E; link+=2*dim_list){
	//		printf("Link: %d, nodes: %d %d \n", link/(2*dim_list), link_list[link], link_list[link+dim_list]);
	//	}
	//	for(int i=0;i<N;i++){
	//		for(int j=0; j<N;j++){
	//			printf("%d ", A[i*N+j]);
	//		}
	//		printf("\n");
	//	}
	//	for(int i=0;i<N;i++){
	//		for(int j=0; j<N;j++){
	//			printf("%1f ", D[i*N+j]);
	//		}
	//		printf("\n");
	//	}
	//	return ;

	int q = 0;
	for (int u=0; u< iterations; u++) {

		int cond = 0;
		while (cond<1) {
			printf("cond: %d, q: %d \n", cond, q);
			q += 1;

			int rand_edge = floor(drand48() * E) ;
			int first_link_index=rand_edge*dim_link_list*2;
			i= link_list[first_link_index];
			j= link_list[first_link_index + dim_link_list];

			//			printf("Length link list: %ld, %d \n", ARR_SIZE(link_list), E);
			//			printf("Link index: %d,  Link i: %d, Link j: %d \n", first_link_index, i ,j );

			float *d_i_all;
			d_i_all=distance(D, i,list_all_neighbors, N, N );
			//			for (int it=0; it<N; it++) {
			//				printf("d_i_all %d - %d: %2f \n ", i, list_all_neighbors[it], d_i_all[it]);
			//			}

			printf("d_i_all\n");
			printFloatArray(d_i_all, N);

			// Create mask and provide distances in mask
			int nb_count=0;
			int mask[N];
			for (int d=0; d<N; d++){
				float Dist_j = d_i_all[d] - d_i_all[j];
				if (fabs(Dist_j) < tolerance * d_i_all[j]){
					mask[d]=1;
					nb_count++;
					//					printf("Node %d, nb_count %d, Dist_j: %2f \n", d, nb_count, Dist_j);
				} else {
					mask[d]=0;
				}
			}

			// We have to exclude the link i<->j from mask! Corresponds to False
			if (nb_count>2){
				mask[i] = 0;
				mask[j] = 0;
			} else {
				continue;
			}

			// Create list of possible neighbors
			int possible_nbs [nb_count-2];
			int tmp=0;
			for (int d=0; d<N; d++){
				if (mask[d]==1){
					possible_nbs[tmp]=d;
					tmp++;
				}
			}

			// Permutation of possible_nbs list
			if (nb_count>2){
				fisher_yates_shuffeling(possible_nbs,nb_count-2);
				printArray(possible_nbs, nb_count-2);
			}

			// Find neighbor in list of possible nbs with correct distance for nb_count -2 times as i, j are not counted)
			l = -1;
			for (int rk=0; rk<nb_count -2 ; rk++) {
				int nk_count=0;
				k=possible_nbs[rk]; // check now the neighbors of k
				int nbs_of_k[N]; // N is the maximum number of neighbor nodes in fully connected network
				for (int link=0; link<2*dim_link_list*E; link+=2*dim_link_list) {
					if (link_list[link]==k ) {
						nbs_of_k[nk_count]=link_list[link+dim_link_list];
						nk_count++;
					}
					// This if is necessary as link list is ordered from smaller node number to larger node number!
					if (link_list[link+dim_link_list]==k ) {
						nbs_of_k[nk_count]=link_list[link];
						nk_count++;
					}
				}

				printf("nbs_of_k k=%d \n", k);
				printArray(nbs_of_k, nk_count);

				// Test if possible neighbor occurs in nbs_of_k as well
				if ((in_array(nbs_of_k, nk_count, i)>-1) | (nk_count == 0) ) {
					//					printf("In Array test link i: %d or nk_count: %d, rk %d , while loop %d, step %d !\n", i, nk_count, rk, q, u);
					continue;
				} else {
					float *d_k_all;
					d_k_all= distance(D,k, nbs_of_k, nk_count,N);
					printf("d_k_all: k %d nk_count %d \n", k, nk_count);

					float *d_j_all;
					d_j_all= distance(D,j, nbs_of_k, nk_count,N);
					printf("d_j_all: j	 %d nk_count %d \n ", j, nk_count);

					int mask2[nk_count];
					int any_candidate=0;
					for(int d=0; d<nk_count;d++) {
						float Dist_k_j=d_k_all[d] - d_j_all[d];
						// This mask is applied furthermore
						// mask2 = np.abs(d_k_all - d_i_all[j]) < tolerance * d_i_all[j]
						float Dist_k_i=d_k_all[d] - d_i_all[j];
						if( (fabs(Dist_k_j)<tolerance*d_k_all[d]) &&( fabs(Dist_k_i) < tolerance*d_i_all[j]) ) {
							mask2[d]=1;
							any_candidate+=1;
						} else {
							mask2[d]=0;
						}
					}

					// Analyze list of possible candidates further
					if (any_candidate > 0) {
						//						printf("Any candidate: %d \n", any_candidate);
						int possible_candidates[any_candidate];
						int candidate_count=0;
						for(int tmp=0; tmp<nk_count; tmp++) {
							if (mask2[tmp]==1) {
								possible_candidates[candidate_count]=nbs_of_k[tmp];
								candidate_count++;
							}
						}

						printf("Possible_candidates: %d \n", candidate_count);
						printArray(possible_candidates, candidate_count);

						// Now check for neighbors of l_candidate
						int l_candidate=possible_candidates[rand() % any_candidate];
						int nl_count=0;
						int nbs_of_l[N];
						for(int link=0; link<2*dim_link_list*E; link+=2*dim_link_list){
							if (link_list[link]==l_candidate)
							{
								nbs_of_l[nl_count]=link_list[link+dim_link_list];
								nl_count++;
							}
							if (link_list[link+dim_link_list]==l_candidate)
							{
								nbs_of_l[nl_count]=link_list[link];
								nl_count++;
							}
						}
						printf("Possible_nl_list: %d \n", nl_count);
						printArray(nbs_of_l, nl_count);

						for (int nl=0; nl<nl_count; nl++) {
							// check if j is not in nbs_of_l
							if (in_array(nbs_of_l, nl_count, j)==-1) {
								l=l_candidate;
								printf("l swap found: l %d j %d \n", l,j);
								rk=nb_count; //  An exchange neighbor is found, leave loop over list of possible neighbors
								break;
							}
						}
					}
					free(d_k_all);
					free(d_j_all);
				}
			}
			free(d_i_all);
			if (l==-1) {
				printf("Return to beginning of while loop! \n");
				continue;  // Returns to beginning of while loop
			} else // Apply change
			{
				printf("Apply changes in adjacency matrix: i %d j %d; k %d l%d \n", i, j, k, l);
				// start changing links in adjacency matrix and link list!
				cond = 1;

				A[i*N + j] =  A[j*N + i] = 0;  // Delete link i<->j
				A[k*N + k] =  A[l*N + k] = 0;  // Delete link k<->l
				A[i*N + k] =  A[k*N + i] = 1;  // Add link i<->k
				A[j*N + l] =  A[l*N + j] = 1;  // Add link j<->l

				printf("Changed Adj matrix. \n");

				// Now find id of second_link_index k<->l
				int second_link_index=0;

				//	for(int link=0; link<2*dim_list*E; link+=2*dim_list){
				//		printf("Link: %d, nodes: %d %d \n", link/(2*dim_list), link_list[link], link_list[link+dim_list]);
				//	}

				for (int link=0; link<2*dim_link_list*E; link+=2*dim_link_list){
					if ( (link_list[link]==k && link_list[link+dim_link_list]==l ) ||  (link_list[link]==l && link_list[link+dim_link_list]==k ) ) {
						second_link_index=link;
						break;
					}
				}
				printf("Second link index found: %d \n", second_link_index);

				// Now update the link list to i<->k j<->l
				link_list[first_link_index] = i;
				link_list[first_link_index+ dim_link_list] = k;
				link_list[second_link_index] = j;
				link_list[second_link_index+dim_link_list] = l;

				printf("link list updated! \n") ;
			}
		}
	}
	printf("Trials %d, Iterations %d \n", q, iterations);
}







void _randomly_rewire_geomodel_I_fast(int iterations, float eps, short *A,
		float *D, int E, int N, int *edges)  {

	int i, s, t, k, l, edge1, edge2, count;
	//int j, neighbor_s_index, neighbor_t_index;
	//int neighbor_k_index, neighbor_l_index;

	//  Create list of neighbors
	//for (int i = 0; i < N; i++) {
	//
	//    count = 0;
	//
	//    for (int j = 0; j < N; j++) {
	//        if (A(i,j) == 1) {
	//            neighbors(i,count) = j;
	//            count++;
	//        }
	//    }
	//}

	//  Initialize random number generator
	srand48(time(0));

	i = 0;
	count = 0;
	while (i < iterations) {
		//  Randomly choose 2 edges
		edge1 = floor(drand48() * E);
		edge2 = floor(drand48() * E);

		s = edges[edge1*N+0];
		t = edges[edge1*N+1];

		k = edges[edge2*N+0];
		l = edges[edge2*N+1];

		//  Randomly choose 2 nodes
		//s = floor(drand48() * N);
		//k = floor(drand48() * N);

		//  Randomly choose 1 neighbor of each
		//neighbor_s_index = floor(drand48() * degree(s));
		//neighbor_k_index = floor(drand48() * degree(k));
		//t = neighbors(s,neighbor_s_index);
		//l = neighbors(k,neighbor_k_index);

		count++;

		//  Proceed only if s != k, s != l, t != k, t != l
		if (s != k && s != l && t != k && t != l) {
			// Proceed only if the new links {s,l} and {t,k}
			// do NOT already exist
			if (A[s*N+l] == 0 && A[t*N+k] == 0) {
				// Proceed only if the link lengths fulfill condition C1
				if ((fabs(D[s*N+t] - D[k*N+t]) < eps &&
						fabs(D[k*N+l] - D[s*N+l]) < eps ) ||
						(fabs(D[s*N+t] - D[s*N+l]) < eps &&
								fabs(D[k*N+l] - D[k*N+t]) < eps )) {
					// Now rewire the links symmetrically
					// and increase i by 1
					A[s*N+t] = A[t*N+s] = 0;
					A[k*N+l] = A[l*N+k] = 0;
					A[s*N+l] = A[l*N+s] = 1;
					A[t*N+k] = A[k*N+t] = 1;

					edges[edge1*N+0] = s;
					edges[edge1*N+1] = l;
					edges[edge2*N+0] = k;
					edges[edge2*N+1] = t;

					//  Update neighbor lists of all 4 involved nodes
					//neighbors(s,neighbor_s_index) = l;
					//neighbors(k,neighbor_k_index) = t;

					//neighbor_t_index = 0;
					//while (neighbors(t,neighbor_t_index) != s) {
					//    neighbor_t_index++;
					//}
					//neighbors(t,neighbor_t_index) = k;

					//neighbor_l_index = 0;
					//while (neighbors(l,neighbor_l_index) != k) {
					//    neighbor_l_index++;
					//}
					//neighbors(l,neighbor_l_index) = s;

					i++;
				}
			}
		}
	}
	printf("Trials %d, Rewirings %d", count, iterations);
}


void _randomly_rewire_geomodel_II_fast(int iterations, float eps, short *A,
		float *D, int E, int N, int *edges)  {

	int i, s, t, k, l, edge1, edge2;

	//  Initialize random number generator
	srand48(time(0));

	i = 0;
	while (i < iterations) {
		//  Randomly choose 2 edges
		edge1 = floor(drand48() * E);
		edge2 = floor(drand48() * E);

		s = edges[edge1*N+0];
		t = edges[edge1*N+1];

		k = edges[edge2*N+0];
		l = edges[edge2*N+1];

		//  Proceed only if s != k, s != l, t != k, t != l
		if (s != k && s != l && t != k && t != l) {
			// Proceed only if the new links {s,l} and {t,k}
			// do NOT already exist
			if (A[s*N+l] == 0 && A[t*N+k] == 0) {

				// Proceed only if the link lengths fulfill condition C2
				if (fabs(D[s*N+t] - D[s*N+l]) < eps &&
						fabs(D[t*N+s] - D[t*N+k]) < eps &&
						fabs(D[k*N+l] - D[k*N+t]) < eps &&
						fabs(D[l*N+k] - D[l*N+s]) < eps ) {
					// Now rewire the links symmetrically
					// and increase i by 1
					A[s*N+t] = A[t*N+s] = 0;
					A[k*N+l] = A[l*N+k] = 0;
					A[s*N+l] = A[l*N+s] = 1;
					A[t*N+k] = A[k*N+t] = 1;

					edges[edge1*N+0] = s;
					edges[edge1*N+1] = l;
					edges[edge2*N+0] = k;
					edges[edge2*N+1] = t;

					i++;
				}
			}
		}
	}
}

void _randomly_rewire_geomodel_III_fast(int iterations, float eps, short *A,
		float *D, int E, int N, int *edges, int *degree)  {

	int i, s, t, k, l, edge1, edge2;

	//  Initialize random number generator
	srand48(time(0));

	i = 0;
	while (i < iterations) {
		//  Randomly choose 2 edges
		edge1 = floor(drand48() * E);
		edge2 = floor(drand48() * E);

		s = edges[edge1*N+0];
		t = edges[edge1*N+1];

		k = edges[edge2*N+0];
		l = edges[edge2*N+1];

		//  Proceed only if s != k, s != l, t != k, t != l
		if (s != k && s != l && t != k && t != l) {
			// Proceed only if the new links {s,l} and {t,k}
			// do NOT already exist
			if (A[s*N+l] == 0 && A[t*N+k] == 0) {
				// Proceed only if degree-degree correlations
				// will not be changed
				if (degree[s] == degree[k] && degree[t] == degree[l]) {
					// Proceed only if the link lengths
					// fulfill condition C2
					if (fabs(D[s*N+t] - D[s*N+l]) < eps &&
							fabs(D[t*N+s] - D[t*N+k]) < eps &&
							fabs(D[k*N+l] - D[k*N+t]) < eps &&
							fabs(D[l*N+k] - D[l*N+s]) < eps ) {
						// Now rewire the links
						// symmetrically and increase i by 1
						A[s*N+t] = A[t*N+s] = 0;
						A[k*N+l] = A[l*N+k] = 0;
						A[s*N+l] = A[l*N+s] = 1;
						A[t*N+k] = A[k*N+t] = 1;

						edges[edge1*N+0] = s;
						edges[edge1*N+1] = l;
						edges[edge2*N+0] = k;
						edges[edge2*N+1] = t;

						i++;
					}
				}
			}
		}
	}
}


double _higher_order_transitivity4_fast(int N, short *A)  {
	long cliques, stars;

	//  Initialize
	cliques = 0;
	stars = 0;
	//  Iterate over all nodes
	for (int v = 0; v < N; v++) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				for (int k = 0; k < N; k++) {
					if (A[v*N+i] == 1 && A[v*N+j] == 1 &&
							A[v*N+k] == 1) {
						stars++;
						if (A[i*N+j] == 1 && A[i*N+k] == 1 &&
								A[j*N+k] == 1)
							cliques++;
					}
				}
			}
		}
	}
	double T = ((double) cliques) / ((double) stars);
	return T;
}


// network ====================================================================

void _do_nsi_hamming_clustering_fast(int n2, int nActiveIndices, float mind0,
		float minwp0, int lastunited, int part1, int part2, float *distances,
		int *theActiveIndices, float *linkedWeights, float *weightProducts,
		float *errors, float *result, int *mayJoin)  {


	int i1, i2, i3, c3;
	int newpart1=0;
	int newpart2=0;
	double d, lw, mind=mind0, minwp=minwp0;
	for (i1=0; i1<nActiveIndices; i1++) {
		int c1 = theActiveIndices[i1];
		if ((lastunited==-1) || (c1==lastunited)) {
			for (i2=0; i2<i1; i2++) {
				int c2 = theActiveIndices[i2];
				if (mayJoin[c1*n2+c2]>0) {
					d = 0.0;
					for (i3=0; i3<i2; i3++) {
						c3 = theActiveIndices[i3];
						lw = linkedWeights[c1*n2+c3]
										   + linkedWeights[c2*n2+c3];
						d += fmin(lw,weightProducts[c1*n2+c3]
													+ weightProducts[c2*n2+c3]-lw)
                                		- errors[c1*n2+c3] - errors[c2*n2+c3];
					}
					for (i3=i2+1; i3<i1; i3++) {
						c3 = theActiveIndices[i3];
						lw = linkedWeights[c1*n2+c3]
										   + linkedWeights[c2*n2+c3];
						d += fmin(lw,weightProducts[c1*n2+c3]
													+ weightProducts[c2*n2+c3]-lw)
                                		- errors[c1*n2+c3] - errors[c2*n2+c3];
					}
					for (i3=i1+1; i3<nActiveIndices; i3++) {
						c3 = theActiveIndices[i3];
						lw = linkedWeights[c1*n2+c3]
										   + linkedWeights[c2*n2+c3];
						d += fmin(lw,weightProducts[c1*n2+c3]
													+ weightProducts[c2*n2+c3]-lw)
                                		- errors[c1*n2+c3] - errors[c2*n2+c3];
					}
					double e = weightProducts[c1*n2+c2]
											  - 2.0*linkedWeights[c1*n2+c2];
					if (e>0.0) d += e;
					distances[c1*n2+c2] = d;
					if ((d<mind) ||
							((d==mind) &&
									(weightProducts[c1*n2+c2]<minwp))) {
						mind = d;
						minwp = weightProducts[c1*n2+c2];
						newpart1 = c1;
						newpart2 = c2;
					}
				}
			}
		} else {
			for (i2=0; i2<i1; i2++) {
				int c2 = theActiveIndices[i2];
				if (mayJoin[c1*n2+c2]>0) {
					double lw_united = linkedWeights[c1*n2+lastunited]
													 + linkedWeights[c2*n2+lastunited],
													 lw_part1 = linkedWeights[c1*n2+part1]
																			  + linkedWeights[c2*n2+part1],
																			  lw_part2 = linkedWeights[c1*n2+part2]
																									   + linkedWeights[c2*n2+part2];
					distances[c1*n2+c2] +=
							(fmin(lw_united, weightProducts[c1*n2+lastunited]
															+ weightProducts[c2*n2+lastunited]
																			 - lw_united)
									- errors[c1*n2+lastunited]
											 - errors[c2*n2+lastunited])
											 - (fmin(lw_part1,weightProducts[c1*n2+part1]
																			 + weightProducts[c2*n2+part1] - lw_part1)
													 - errors[c1*n2+part1] - errors[c2*n2+part1])
													 - (fmin(lw_part2,weightProducts[c1*n2+part2]
																					 + weightProducts[c2*n2+part2] -lw_part2)
															 - errors[c1*n2+part2] - errors[c2*n2+part2]);
					d = distances[c1*n2+c2];
					if ((d<mind) ||
							((d==mind) &&
									(weightProducts[c1*n2+c2]<minwp))) {
						mind = d;
						minwp = weightProducts[c1*n2+c2];
						newpart1 = c1;
						newpart2 = c2;
					}
				}
			}
		}
	}
	result[0] = mind;
	result[1] = newpart1;
	result[2] = newpart2;
}


double _vertex_current_flow_betweenness_fast(int N, double Is, double It,
		float *admittance, float *R, int i) {

	double VCFB=0.0;
	int t=0;
	int s=0;
	int j=0;
	double I=0;

	for(t=0;t<N;t++){
		for(s=0; s<t; s++){
			I = 0.0;
			if(i == t || i == s){
				continue;
			}
			else{
				for(j=0;j<N;j++){
					I += admittance[i*N+j]*
							fabs( Is*(R[i*N+s]-R[j*N+s])+
									It*(R[j*N+t]-R[i*N+t])
							) / 2.0;
				} // for  j
			}
			VCFB += 2.0*I/(N*(N-1));
		} // for s
	} // for t

	return VCFB;
}


void _edge_current_flow_betweenness_fast(int N, double Is, double It, 
		float *admittance, float *R, float *ECFB) {

	int i=0;
	int j=0;
	int t=0;
	int s=0;
	double I = 0.0;

	for(i=0; i<N; i++){
		for(j=0;j<N;j++){
			I = 0.0;
			for(t=0;t<N;t++){
				for(s=0; s<t; s++){
					I += admittance[i*N+j]*\
							fabs(Is*(R[i*N+s]-R[j*N+s])+
									It*(R[j*N+t]-R[i*N+t]));
				} //for s
			} // for t
			ECFB[i*N+j] += 2.*I/(N*(N-1));
		} // for j
	} // for i
}
