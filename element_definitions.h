/*

Copyright (C) 2003 The GeoFramework Consortium

This file is part of Ellipsis3D.

Ellipsis3D is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License, version 2,
as published by the Free Software Foundation.

Ellipsis3D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Authors:
  Louis Moresi <louis.moresi@sci.monash.edu>
  Richard Albert <richard.albert@exxonmobil.com>

*/



/*	Header file containing the appropriate definitions for the elements being used 
in the problem.   Include in all of the files linking in to the element program.  */

/* #define const  */

#define ELNMAX 27
#define ELN 8  /* for this element, the max nodes/vpoints/ppoints, used in indexing */
#define ELV 8
#define ELP 1
#define EL1N 4   /*RAA: 10/5/01, ELN and ELV are 8 for 3D, while EL1N and EL1V should be 4 */
#define EL1V 4
#define EL1P 1
#define GNVI ELV*ELN
#define GNVI2 ELV*ELN*2
#define GNVI3 ELV*ELN*3
#define GNVI6 ELV*ELN*6
#define GNPI ELP*ELN
#define GNPI2 ELP*ELN*2
#define GNPI3 ELP*ELN*3
#define GNPI6 ELP*ELN*6
#define GN1VI EL1V*EL1N
#define GN1PI EL1P*EL1N

#define GNVINDEX(n,v) ((ELV*((n)-1)) + ((v)-1))
#define GNPINDEX(n,p) ((ELP*((n)-1)) + ((p)-1))

#define GNVXINDEX(d,n,v) ((ELV*ELN*(d))+(ELV*((n)-1)) + ((v)-1))
#define GNPXINDEX(d,n,p) ((ELP*ELN*(d))+(ELP*((n)-1)) + ((p)-1))
#define GNVXSHORT(d,i) ((ELV*ELN*(d))+(i))
#define GNPXSHORT(d,i) ((ELP*ELN*(d))+(i))

#define TWENTY_SEVEN_NODES_CUBIC 27
#define TWENTY_SIX_NODES_CUBIC 26
#define TWENTY_NODES_CUBIC 20
#define EIGHT_NODES_CUBIC 8
#define NINE_NODES_QUAD 9
#define EIGHT_NODES_QUAD 8
#define SEVEN_NODES_QUAD 7
#define SIX_NODES_QUAD 6
#define FIVE_NODES_QUAD 5
#define FOUR_NODES_QUAD 4
#define THREE_NODES_QUAD 3
#define ONE_NODE_QUAD 1
#define ONE_NODE_CUBIC 1

/*
#define GNVINDEX(n,v) ((ELN*((v)-1)) + (n-1))
#define GNPINDEX(n,p) ((ELN*((p)-1)) + (n-1))
#define GNVXINDEX(d,n,v) ((ELV*ELN*(d))+(ELN*((v)-1)) + (n-1))
#define GNPXINDEX(d,n,p) ((ELP*ELN*(d))+(ELN*((p)-1)) + (n-1))
*/


#define GMVINDEX(n,v) ((EL1V*((n)-1)) + (v-1))
#define GMPINDEX(n,p) ((EL1P*((n)-1)) + (p-1))
#define GMVXINDEX(d,n,v) ((EL1V*EL1N*(d))+(EL1V*((n)-1)) + (v-1))
#define GMPXINDEX(d,n,p) ((EL1P*EL1N*(d))+(EL1P*((n)-1)) + (p-1))

#define GMVGAMMA(i,n) (4*(i) + (n))


#define BVINDEX(di,dj,n,v) ((GNVI3*(di-1))+(GNVI*(dj-1))+(ELV*((n)-1)) + (v-1))
#define BVXINDEX(di,dj,i,n,v) ((GNVI6*(di-1))+(GNVI2*(dj-1))+(GNVI*(i-1))+(ELV*((n)-1)) + (v-1))

#define BPINDEX(di,dj,n,v) ((GNPI3*(di-1))+(GNPI*(dj-1))+(ELP*((n)-1)) + (v-1))
#define BPXINDEX(di,dj,i,n,v) ((GNPI6*(di-1))+(GNPI2*(dj-1))+(GNPI*(i-1))+(ELP*((n)-1)) + (v-1))

/* Element definitions */

static const int enodes[4] = {0,2,4,8};
static const int pnodes[4] = {0,1,1,1};
static const int vpoints[4] = {0,2,4,8};
static const int ppoints[4] = {0,1,1,1};
static const int onedvpoints[4] = {0,0,2,4};
static const int onedppoints[4] = {0,0,1,1};

/* More cumbersome, but more optimizable versions ! */

#define ENODES3D 8
#define ENODES2D 4
#define VPOINTS3D 8
#define VPOINTS2D 4
#define PPOINTS3D 1
#define PPOINTS2D 1
#define V1DPOINTS3D 4
#define V1DPOINTS2D 2
#define P1DPOINTS3D 1
#define P1DPOINTS2D 1

/*  As we use arrays starting at index 1  */
static const int sfnarraysize[] = {0,3,5,9};
static const int onedsfnarsize[] = {0,0,3,5};
static const int gptarraysize[] = {0,3,5,9};
static const int pptarraysize[] = {0,2,2,2};
static const int node_array_size[] = {0,3,5,9};
static const int pnode_array_size[] = {0,2,2,2};
static const int onedgptarsize[] = {0,0,3,5};
static const int onedpptarsize[] = {0,0,2,2};
/*  As we use arrays starting at index 0  */
static const int loc_mat_size[4][3] = {{0,0,0},{0,4,4},{0,8,12},{0,24,48}} ;
static const int stored_mat_size[] = {0,10,36,300};
static const int node_mat_size[] = {0,3,18,81};
static const int ploc_mat_size[] = {0,1,1,1};
/*  Inter-node information */
static const int max_node_interaction[] = {0,3,9,27};
static const int max_eqn_interaction[]={0,3,18,81};
static const int max_els_per_node[]={0,2,4,8};

#define MAX_EQN_INT_3D 81
#define MAX_EQN_INT_2D 18
#define MAX_EQN_INT_2pt5D 27
#define MAX_NODE_INT_3D 27
#define MAX_NODE_INT_2D 9
#define MAX_NODE_INT_2pt5D 9

static const int seq[]={0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153,171,190,210,231,253,276,300};


#define BETA  0.57735026918962576451 /* (sqrt3)/3 */


/*  4-pt co-ords  */

/*   Declare everything to be static: no changes are made to this data ... it is 
easiest to redeclare it all each time and have just one file of definitions   */

static const struct One_d_int_points {
	float x[2];
	float weight[3]; 	 } 
			g_1d[5] = 
                        {	{{0.0,0.0},{0.0,0.0,0.0}},
				{{-BETA,-BETA},{0.0,1.0,1.0}},
				{{-BETA, BETA},{0.0,1.0,1.0}}, 
				{{ BETA, BETA},{0.0,1.0,1.0}}, 
				{{ BETA,-BETA},{0.0,1.0,1.0}}   } ;

static const struct One_d_int_points 
			p_1d[2] = 
			{	{{0.0,0.0},{0.0,0.0,0.0}},
				{{0.0,0.0},{0.0,2.0,4.0}} };
	                                                  	   
static const struct Int_points {
	float x[3]; 
	float weight[3];    }
                        g_point[9] = 
			{	{{0.0,0.0,0.0},{0.0,0.0,0.0}},
				{{-BETA,-BETA,-BETA},{0.0,1.0,1.0}},
				{{-BETA, BETA,-BETA},{0.0,1.0,1.0}},
                                {{ BETA, BETA,-BETA},{0.0,1.0,1.0}},
                                {{ BETA,-BETA,-BETA},{0.0,1.0,1.0}},
                                {{-BETA,-BETA, BETA},{0.0,1.0,1.0}},
				{{-BETA, BETA, BETA},{0.0,1.0,1.0}},
                                {{ BETA, BETA, BETA},{0.0,1.0,1.0}},
                                {{ BETA,-BETA, BETA},{0.0,1.0,1.0}} };
  
  
static const struct Int_points 
			s_point[5] = 
			{	{{0.0,0.0,0.0},{0.0,0.0,0.0}},
				{{-BETA,0.0,-BETA},{0.0,2.0,2.0}},
				{{ BETA,0.0,-BETA},{0.0,2.0,2.0}},
                                {{ BETA,0.0, BETA},{0.0,2.0,2.0}},
                                {{-BETA,0.0, BETA},{0.0,2.0,2.0}} };
  
  
static const struct Int_points 
			p_point[2] =
			{	{{0.0,0.0,0.0},{0.0,0.0,0.0}},
				{{0.0,0.0,0.0},{0.0,4.0,8.0}} };
	                  
static const struct Node_points { 
	float x[3]; }
			node_point[9] =
			{	{0.0,0.0,0.0},
				{-1.0,-1.0,-1.0},
				{-1.0, 1.0,-1.0},
				{ 1.0, 1.0,-1.0},
				{ 1.0,-1.0,-1.0},
				{-1.0,-1.0, 1.0},
				{-1.0, 1.0, 1.0},
				{ 1.0, 1.0, 1.0},
				{ 1.0,-1.0, 1.0}	};
				
static const struct Internal_structure  /* integer coordinates relative to node 1 */
	{	int vector[3]; } offset[9] =
				{	{0,0,0},
					{0,0,0},
					{0,1,0},
					{1,1,0},
					{1,0,0},
					{0,0,1},
					{0,1,1},
					{1,1,1},
					{1,0,1}	};
	
static const struct E_loc_info {
	int minus[3];
	int plus[3];
	float delta[3];
	int num_nebrs[3];
	int node_nebrs[3][2]; } 
		loc[9] = {
			{{0,0,0},{0,0,0},{0.0,0.0,0.0},{0,0,0},{{0,0},{0,0},{0,0}}},
			{{0,0,0},{4,2,5},{0.5,0.5,0.5},{2,2,2},{{1,4},{1,2},{1,5}}},/* 1 */
			{{0,1,0},{3,0,6},{0.5,0.5,0.5},{2,2,2},{{2,3},{1,2},{2,6}}},
			{{2,4,0},{0,0,7},{0.5,0.5,0.5},{2,2,2},{{2,3},{4,3},{3,7}}},
			{{1,0,0},{0,3,8},{0.5,0.5,0.5},{2,2,2},{{1,4},{4,3},{4,8}}},
	 		{{0,0,1},{8,6,0},{0.5,0.5,0.5},{2,2,2},{{5,8},{5,6},{1,5}}}, /* 5 */
			{{0,5,2},{7,0,0},{0.5,0.5,0.5},{2,2,2},{{6,7},{5,6},{2,6}}},
			{{6,8,3},{0,0,0},{0.5,0.5,0.5},{2,2,2},{{6,7},{8,7},{3,7}}},
			{{5,0,4},{0,7,0},{0.5,0.5,0.5},{2,2,2},{{5,8},{8,7},{4,8}}} };
			
	
static const int bb[3][9] = {	{0,1,1,2,2,1,1,2,2},  /* x dirn */
			        {0,1,2,2,1,1,2,2,1},  /* z dirn */
			        {0,1,1,1,1,2,2,2,2}   /* y dirn */ };

static const int ccc[3][5] = { {0,1,3,0,0},    /* use these nodes 1d */
		               {0,1,2,3,4},
		               {0,0,0,0,0}}; /* ok ordering 2d */

static const int ee[3][2] = { {0,1},	/* x dirn */
			{0,1},
			{0,1} };


/* I am node %d in an element of the local group, which element is it .... */

static const int loc_e_stencil[3][9] = {
    {0,2,1,0,0,0,0,0,0},
    {0,4,3,1,2,0,0,0,0},
    {0,8,7,5,6,4,3,1,2} } ;
