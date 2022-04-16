import math
import random
from copy import deepcopy
import matplotlib.pyplot as plt


class Chromosome:
    __repr_length = 0
    __interval = (0, 0)
    __precision = 0

    def __init__(self, interval, precision):
        self.__interval = interval
        self.__precision = precision
        if Chromosome.__repr_length == 0:
            Chromosome.__repr_length = Chromosome.calculate_length(interval, precision)
        self.genes = [random.randint(0, 1) for _ in range(Chromosome.__repr_length)]

    @staticmethod
    def calculate_length(interval, p):
        a, b = interval
        power_of_10 = 10 ** p
        return int(math.log2((b - a) * power_of_10)) + 1

    def convert_base10(self):
        power = 1
        result = 0
        for value in reversed(self.genes):
            result += power * value
            power *= 2
        result = (self.__interval[1] - self.__interval[0]) / (2 ** Chromosome.__repr_length - 1) * \
                 result + self.__interval[0]

        return round(result, self.__precision)  # posibil sa trebuiasca schimbat

    def set_genes(self, genes):
        self.genes = deepcopy(genes)

    @staticmethod
    def get_length():
        return Chromosome.__repr_length

    def __len__(self):
        return len(self.genes)

    def __setitem__(self, key, value):
        self.genes[key] = value

    def __getitem__(self, item):
        return self.genes[item]

    def __str__(self):
        return "".join([str(gene) for gene in self.genes])


class Population:
    """
        size - the number of individuals in the Population
        interval - the definition domain of the function
        crossover - the probability that two chromosomes will combine
        mutation - the probability that a chromosome will mutate
        params - the parameters a, b, c of the function
                    f = a * x ^ 2 + b * x + c we want to maximise
        subjects - a vector of chromosomes representing our sample
    """
    __size = 0
    __interval = (0, 0)
    __precision = 0
    __crossover = 0
    __mutation = 0
    __params = (0, 0, 0)

    def __init__(self, size, interval, params, precision, crossover, mutation):
        Population.__size = size
        Population.__interval = interval
        Population.__params = params
        Population.__precision = precision
        Population.__crossover = crossover
        Population.__mutation = mutation
        self.subjects = [Chromosome(interval, precision) for _ in range(size)]
        self.subjects = sorted(self.subjects, key=lambda x: self.calculate_f(x))
        self.printing = True
        self.max_value = [Population.calculate_f(self.subjects[Population.__size - 1])]

    def set_subjects(self, subjects):
        obj = deepcopy(subjects)
        self.subjects = obj

    def get_subject_by_index(self, idx):
        return self.subjects[idx]

    @staticmethod
    def calculate_f(chromosome):
        x = chromosome.convert_base10()
        a, b, c = Population.__params
        f = a * (x ** 2) + b * x + c
        return round(f, Population.__precision)

    def calculate_probability_of_selection(self):
        f_values = [self.calculate_f(subject) for subject in self.subjects]
        s = sum(f_values)
        return [round(f_value / s, Population.__precision) for f_value in f_values]

    def generate_selection_intervals(self):
        ps = self.calculate_probability_of_selection()

        if self.printing:
            f.write("Probabilities of selection\n")
            for (i, p) in enumerate(ps):
                f.write(f"Chromosome {i + 1}    probability: {p}\n")
            f.write("\n")

        left = 0
        right = ps[0]
        intervals = []
        for value in ps[1:]:
            intervals.append((left, right))
            left = right
            right = round((right + value), Population.__precision)
        right = 1
        intervals.append((left, right))
        return intervals

    def find_interval(self, left, right, target, intervals):
        if target < intervals[left][1]:
            return left
        if target >= intervals[right][0]:
            return right

        m = (right + left) // 2
        interval = intervals[m]

        if interval[0] <= target < interval[1]:
            return m

        if target < interval[0]:
            return self.find_interval(left, m - 1, target, intervals)
        else:
            return self.find_interval(m + 1, right, target, intervals)

    @staticmethod
    def merge_chromosomes(subjects, printing):
        shuffled_sub = [elem for elem in enumerate(subjects)]
        random.shuffle(shuffled_sub)
        new_subjects = []
        if printing:
            f.write("\n")
        # if the number of chromosomes is odd, the last chromosome remains unchanged
        for i in range(0, len(subjects) - 1, 2):
            idx1, chromosome1 = shuffled_sub[i]
            idx2, chromosome2 = shuffled_sub[i + 1]

            poz = random.randint(1, Chromosome.get_length() - 2)
            new_genes1 = chromosome1[:poz] + chromosome2[poz:]
            new_chromosome1 = Chromosome(Population.__interval, Population.__precision)
            new_chromosome1.set_genes(new_genes1)
            new_subjects.append(new_chromosome1)

            new_genes2 = chromosome2[:poz] + chromosome1[poz:]
            new_chromosome2 = Chromosome(Population.__interval, Population.__precision)
            new_chromosome2.set_genes(new_genes2)
            new_subjects.append(new_chromosome2)

            if printing:
                f.write(f"Crossing over between chromosome {idx1} and {idx2} with breaking point {poz}:\n")
                f.write(str(new_chromosome1) + "  " + str(new_chromosome2) + "\n")

        if len(subjects) % 2 == 1:
            new_subjects.append(subjects[len(subjects) - 1])

        return new_subjects

    @staticmethod
    def crossing_over(subjects, printing=False):
        will_merge = []
        selected_indexes = set()
        if printing:
            f.write(f"Probability of crossing-over {Population.__crossover}\n")
        for (idx, chromosome) in enumerate(subjects):
            num = round(random.random(), Population.__precision)
            if num < Population.__crossover:
                will_merge.append(chromosome)
                selected_indexes.add(idx)
            if printing:
                f.write(f"{idx + 1}: {chromosome}    u= {num}")
                if num < Population.__crossover:
                    f.write(f" <{Population.__crossover} participates\n")
                else:
                    f.write("\n")

        return Population.merge_chromosomes(will_merge, printing) + \
               [subjects[i] for i in range(len(subjects)) if i not in selected_indexes]

    @staticmethod
    def mutate(subjects, printing):
        num = int(Population.__size * Population.__mutation)
        idxs = set()

        for i in range(num):
            idx = random.randint(0, len(subjects) - 1)
            idxs.add(idx)
            for j in range(len(subjects[idx])):
                if random.random() > Population.__mutation:
                    subjects[idx][j] = 1 - subjects[idx][j]
        if printing:
            f.write(f"Probability of mutation {Population.__mutation}\n")
            f.write(f"Modified chromosomes: \n")
            for idx in idxs:
                f.write(str(idx) + "\n")
        return subjects

    def selection(self):

        if self.printing:
            f.write("Initial population:\n")
            for i in range(self.__size):
                chromosome = self.subjects[i]
                f.write(f"\t {i + 1}: {chromosome} x= {chromosome.convert_base10()} "
                        f"f= {self.calculate_f(chromosome)}\n")
            f.write("\n")

        intervals = self.generate_selection_intervals()

        if self.printing:
            f.write("Intervals of probability of selection:\n")
            for interval in intervals:
                f.write(str(interval[0]) + " ")
            f.write("1.0\n")

        # the chromosome with the best fitness goes straight to the next generation
        next_generation = [deepcopy(self.subjects[len(self.subjects) - 1])]
        # choose the chromosome using the intervals
        first_phase = []
        for i in range(Population.__size - 1):
            number = round(random.random(), Population.__precision)
            idx = self.find_interval(0, len(intervals) - 1, number, intervals)

            if self.printing:
                f.write(f"u= {number}    select chromosome {idx + 1}\n")

            first_phase.append(self.get_subject_by_index(idx))

        if self.printing:
            f.write("\nAfter selection:\n")
            for i in range(Population.__size - 1):
                chromosome = first_phase[i]
                f.write(f"\t {i + 1}: {chromosome} x= {chromosome.convert_base10()} "
                        f"f= {self.calculate_f(chromosome)}\n")
            f.write(f"\t {Population.__size}: {next_generation[0]} x= {next_generation[0].convert_base10()} "
                        f"f= {self.calculate_f(next_generation[0])}\n")
            f.write("\n")

        # crossing over
        second_phase = Population.crossing_over(first_phase, self.printing)
        if self.printing:
            f.write("\nAfter recombination:\n")
            for i in range(len(second_phase)):
                chromosome = second_phase[i]
                f.write(f"\t {i + 1}: {chromosome} x= {chromosome.convert_base10()} "
                        f"f= {self.calculate_f(chromosome)}\n")
            f.write(f"\t {Population.__size}: {next_generation[0]} x= {next_generation[0].convert_base10()} "
                        f"f= {self.calculate_f(next_generation[0])}\n")
            f.write("\n")
        # mutation
        third_phase = Population.mutate(second_phase, self.printing)
        next_generation.extend(third_phase)
        next_generation = sorted(next_generation, key=lambda x: self.calculate_f(x))
        if self.printing:
            f.write("\nAfter mutation:\n")
            for i in range(len(next_generation)):
                chromosome = next_generation[i]
                f.write(f"\t {i + 1}: {chromosome} x= {chromosome.convert_base10()} "
                        f"f= {self.calculate_f(chromosome)}\n")
            f.write("\nEvolution of the maximum value:\n")
        f.write(str(Population.calculate_f(next_generation[Population.__size - 1])) + '\n')
        # update generation
        self.printing = False
        self.max_value.append(Population.calculate_f(next_generation[Population.__size - 1]))
        self.set_subjects(next_generation)

    def __repr__(self):
        s = ""
        for (idx, individual) in enumerate(self.subjects):
            line = str(idx) + ": " + str(individual) + \
                   " x= " + str(individual.convert_base10()) + \
                   " f= " + str(Population.calculate_f(individual)) + '\n'
            s += line
        return s


def parse_input(filename):
    f = open(filename)
    size = int(f.readline())
    interval = tuple(float(x) for x in f.readline().split())
    params = tuple(float(x) for x in f.readline().split())
    precision = int(f.readline())
    crossover = float(f.readline())
    mutation = float(f.readline())
    steps = int(f.readline())
    return (steps, Population(size, interval, params, precision, crossover, mutation))


(steps, population) = parse_input("input.txt")
f = open("output.txt", "w")

for i in range(steps):
    population.selection()
f.write(str(Population.calculate_f(population.subjects[len(population.subjects) - 1])) + '\n')
f.close()

plt.plot(population.max_value)
plt.draw()
plt.show()
