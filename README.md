# genome
This is a c++ library that allows you to use genetic algorithms very easily.Now create an ai for your games with the power of evolution at your hand, introducing genome.h

Documentation for genome
=========================

Author
------
Abhinav - p (@AIFactory)

Introduction
------------

genome is a genetic algorithm library made for c++. It can be used to automate games and making evolutionary simulations using a c++ program. It uses the power of evolution and neural network to reach an optimal solution to a given problem provided sufficient data about the problem is supplied.

Requirements
------------

genome requires the latest stable version of c++ installed in your computer.

genome is a library that works on top of c++ to simulate genetic evolution.

NOTE :(SFML) genome includes a function to visualise the neural network graphically. Since the graphics used for rendering the neural network is based on a library called SFML.

SFML has to be installed and properly working in your computer inorder to use those functions

NOTE: Neural network visualisation functions are commented out by default and are inaccessible. Uncomment the specified section of code to use visualisation functions only if your computer has SFML installed and working properly.


How to include
-------------

#include "genome.h"


Getting started with genome library
------------------------------------

genome libary has a class called genetic where all the functions for genetic algorithms are present.

Create an object of the genetic class to use genome library.

#include "genome.h"

genetic evo(int no_of_inputs,int no_of_outputs,float mutation_rate=0.01,long int population_size=100);

parameters:

evo is an object of class genome

no_of_inputs -> it specifies the number of inputs supplied to the genome library .

no_of_outputs -> it specified the number of outputs to be taken from the genome library.

mutation_rate -> it specifies the mutation rate by which genetic algorithm takes place.

population_size -> size of the population in one generation.

note: 	mutation rate ?-> In genetic algorithm the dna of some species are randomly changed to simulate mutation in the real world. mutation rate is a crucial factor which determine how many species have to be mutated in a given population.

A mutation rate of 0.01 means there is a 1% chance that a species will be mutated after the evolution process.

A higher mutation rate destroys the evolution there by increasing randomness.
A lower mutation rate prevents the species from exploring new areas thus preventing improvement of evolution.

public elements in genetic class
--------------------------------

1)population-> an array of objects having genetic data about the population of the species

public elements accessible from population
------------------------------------------

1)mutation_rate
---------------
evo.population[1].mutation_rate;
this represents the mutation rate of the 2nd agent of the population

2)fitness
---------
evo.population[1].fitness=10;
this refers to the fitness of a species-2 of the population. 
we can set the population for each species like this


Fuctions under genetic class of genome
--------------------------------------

1)setfitness(vector<int> fit_arr)
---------------------------------

=>parameters:

fit_arr -> an array of fitness values for all the agents of the population in the order of the agents in the popluation
note : the length of the fit_arr must match the population_size

=>functionality:
sets the fitness of the species of a generation

=>return type:
void


2)set_population_size(int size)
-------------------------------

=>parameters:

size-> an integer representing population size.

=>functionality:

can be used to change the population size even when the algorithm starts running. changing population size after a specific generation, conditional changing of population etc.

=>return type:
void

3)set_mutation_rate(float mutation)
------------------------------------

=>parameters:

mutation-> mutation rate of the genetic evolution. a value between 0 and 1

=>functionality:

can be used to set and change the mutation rate in the run time of the program. can be used to deacy the mutation_rate in the genetic evolution

=>return type:

void

4)decay_mutation(float decay_rate)
----------------------------------
=>parameters:

decay_rate-> it specifies the decay constant of the mutation, usually nearer to 1 example (0.999) 

=>functionality:

used to decay the mutation rate over generation.-> it means to reduce the mutation rate as the generation increases. this readuces the randomness as the generation increases and thus avoid random actions in a long run of the  program

=>return type:

void

5)get_weights(int agent)
-------------------------

=>parameters:

agent-> an integer representing the agent (0->population_size-1)

=>functionality:

return an array of weights of the agent

=>return type:

vector<float>

6)best_fit_agent()
------------------

=>parameters:
none

=>functionality:

returns the index of the best fit agent of the generation

=>return type:

int


7)get_genes(int index)
-----------------------

=>parameters:

index-> the index we want to get the genes

=>functionality:

returns the array of genes of the agent of the specified index of the current generation.

=>return type:

vector<float>

8)output(vector<double> input, int agent)
-----------------------------------------


=>parameters:

input-> an array of type double having the input values to the algorithm, the input array must have same size as the number of inputs suplied earlier.

agent-> it is agent which we are supplying the input

=>example:
If we are providing the input for 2nd agent and our inputs are x,y,z
output({x,y,z},3);

=>functionality:
returns an output from the genetic algorithm .
return is an interger values starting from 0,
for an output_no of 3 the output function can return 0,1,or 2 depending on the neural network.

=>return:

int 
representing each output

=>code sample:

int action=output({x,y,z},3);
if(action==0){
	jump();
}
else if(action==1){
	duck();
}


9)start_mating()
----------------

=>parameters:
none

=>functionality:

starts the mating process or starts the genetic aglorithm , use it only after setting fitness of all the species of the first generation.
usage without seting a base fitness can lead to errors.

=>return:
int-> it returns the  current generation number being processed

=>code sample:
int gen=start_mating();


10)show_network()
-------------------

=>parameter:
none

=>functionality:
shows a tabluar repesentation of the neural network and the information about the current generation
->neural network of the best fit agent of the generation is shown

=>return:
void

11)get_generation()
-------------------

=>parameter:
none

=>functionality:
returns the current generation that is being processed

=>return:
int

12)all_weights()
------------------

=>parameter:
none

=>functionality:
returns a 2 diamensional array fo weights used to evaluate the output of the neural network.
the return has a size of =(output_no)*(input_no+1)

=>return:

vector<vector<float>>

13)get_inputs()
----------------

=>functionality:

returns the current inputs supplied to the neural network. Can be used as a debugging tool.returns the array of inputs.

=>return:	
vector<double>


global functions in the genome library
---------------------------------------

1)range_map(double var, double xi, double xf, double Xi, double Xf)
-------------------------------------------------------------------

=>parameters:

var-> the variable to be mapped
xi-> the lowest value of variable scale
xf-> the largest value of variable scale
Xi-> the lowest value of the new scale 
Xf-> the largest value of the new scale

=>functionality:

it takes in a value between a range and gives a value scaled between another range.
returns the new_value in new scale

=>return:
double


Functions for neural network visualisation
------------------------------------------

NOTE: inorder to use these functions you need to uncomment a section of the code in the genome.h header file
these functions require the sfml library installed in your system.

1)draw_network(int output_no, int input_no, vector<vector<float>> out_w, vector<float> positions = { 400,100,500,130,50,50 })

----------------------------------------------------------------------

=>parameters:

output_no -> number of outputs
input_no -> number of inputs

out_w -> weights of the best fit agent of generation . Use best_fit_agent() function  here

positions-> array of positions of the neural network to be displayed on screen
			{input_sec_x,input_sec_y,output_sec_x,output_sec_y,input_gap,output_gap}
			sec->section
			x->x cordinate
			y->y cordinate


=>functionality:

return an array of sf::CircleShape objects for drawing circles nodes of the network in the main function

=>return type:

vector<sf::CircleShape>

=>example draw in main function:

vector<sf::CircleShape> net_cir;
net_cir = draw_network(2, 3, g.all_weights(), { 300,100,600,160,100,100 });
for (int c = 0; c < net_cir.size(); c++) {
    window.draw(net_cir[c]);
}


2)draw_network_lines(int output_no, int input_no, vector<vector<float>> out_w, vector<float> positions = { 400,100,500,130,50,50 })

-----------------------------------------------------------------------

=>parameters:

output_no -> number of outputs
input_no -> number of inputs

out_w -> weights of the best fit agent of generation . Use best_fit_agent() function  here

positions-> array of positions of the neural network to be displayed on screen
			{input_sec_x,input_sec_y,output_sec_x,output_sec_y,input_gap,output_gap}
			sec->section
			x->x cordinate
			y->y cordinate


=>functionality:

return an array of sf::RectangleShape objects for drawing line in the main function

=>return type:

vector<sf::RectangleShape>

=>example draw in main function:

vector<sf::RectangleShape> net_line;
net_line = draw_network_lines(2, 3, g.all_weights(), { 300,100,600,160,100,100 });
for (int c = 0; c < net_line.size(); c++) {

    window.draw(net_line[c]);
}
window.display();


about this library
------------------

This library requires no additional libraries except for the neural network visual function.
neural network visual functions works on SFML which has its own standard liscense.

This library is created by Abhinav-P (@ai-factoy) 
licensed under [MIT License](LICENSE)

comment if you are facing any issue with this library.

Hope this library helps .
