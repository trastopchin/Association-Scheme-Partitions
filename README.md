# Association-Scheme-Partitions

This project uses GAP to compute equitable partitions of association schemes. This project uses and extends Izumi Miyamoto and Akihide Hanaki's code that already answers questions about association schemes to help compute equitable partitions. Created as a final research project for our Algebraic Graph Theory course at Grinnell College the spring of 2020 with professor Christopher French.

This project is comprised of the code that we wrote and the paper we wrote about the project. All of the code that we wrote is in the [equitable_partitions.gap](https://github.com/trastopchin/Association-Scheme-Partitions/blob/master/equitable_partitions.gap) file, where we wrote helper methods and algorithms to compute and log equitable partitions of association schemes. The [paper](insert link) we wrote gives some mathematical background on association schemes and equitable partitions, details our programming conventions, documents our methods, and explains our algorithms.

## How to use

To use this project, you must first install the latest version of [GAP](https://www.gap-system.org/). Once you clone this repository, all you need to do is start a GAP session in that directory and read in the file with the command `Read("equitable_partitions.gap");` Then, you have access to each of the methods and algorithms that we have written. For more extensive documentation and explanation about each of the algorithms (other than the in-line documentation), consult the methods section of the [paper](link-to-paper) we wrote about this project.

## Implementation Details

* The methods written (IsIsomorphic, IsEquitablePartition, EquitablePartitions and EquitablePartitionsFast) all rely on exhaustive searches.

* EquitablePartitionsFast is implemented using a hash table in an attempt to speed it up and
not rely on the exhaustive IsIsomorphic and IsEquitablePartition as much.

* More details can be found in the paper linked above.

## Built With

* [GAP](https://www.gap-system.org/) - GAP - Groups, Algorithms, Programming -
a System for Computational Discrete Algebra
* [Classification of Association Schemes](http://math.shinshu-u.ac.jp/~hanaki/as/) - Miyamoto and Hanaki's classification of association schemes

## Authors

* **Gabby Masini** - [GitHub](https://github.com/masiniga)
* **Tal Rastopchin** - [GitHub](https://github.com/trastopchin)
