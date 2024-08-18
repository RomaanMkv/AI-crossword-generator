import copy
import random
from typing import List, Dict, Tuple
import os

# number of inputs files - number of files in inputs directory
NUMBER_OF_CROSSWORDS = len(os.listdir("inputs"))

class Word:
    """
    This class represents the word
    """

    def __init__(self, name):
        # the word itself
        self.name = name
        # horizontal (along OX axis) ir vertical (along OY axis)
        self.position = -1
        # represents the words which are crossed by this word: 0 element - index of the letter of this word, 1 element -
        # index of another word, 2 element - index of the letter in another word
        self.neighbours = []
        # this is used for DFS algorithm
        self.visited = False
        # coordinates first letter of the word in the grid
        self.coordinates = [-1, -1]


class Sample:
    """
    This class represents crossword
    """

    def __init__(self, word_list, dictionary, grid_list, fitness_func):
        # list of all words in crossword
        self.word_list = copy.deepcopy(word_list)
        # dictionary of the crossword
        self.dictionary = copy.deepcopy(dictionary)
        # list of all grids of the crossword (list is needed for different component connectivity)
        self.grid_list = copy.deepcopy(grid_list)
        # fitness function value of this crossword
        self.fitness_function_value = fitness_func


def identical_elements(arr):
    """
    This function determine if all elements of the list are identical
    :param arr: given list
    :return: True if all elements of array are identical, otherwise - False
    """
    for i in range(len(arr) - 1):
        if arr[i] != arr[i + 1]:
            return False
    return True


def get_list_of_component(index, crossword_words: List[Word], lst: List[int]):
    """
    This function determines indexes of all words in connectivity component
    This function uses DFS algorithm
    :param index: index of the word
    :param crossword_words: list of all words of the crossword
    :param lst: list of indexes of this component
    :return: list of indexes of all words in the connectivity component
    """
    if crossword_words[index].visited:
        return lst
    else:
        crossword_words[index].visited = True
        lst.append(index)

    # recursively call this function for all his 'neighbours'
    for nb in crossword_words[index].neighbours:
        lst = get_list_of_component(nb[1], crossword_words, lst)

    return lst


def get_list_of_all_components(crossword_words: List[Word]) -> List[List[int]]:
    """
    This function determines list of all components of the crossword
    :param crossword_words: list of all crossword words
    :return: List of each component words' indexes
    """
    global NUMBER_OF_WORDS

    component_list = []
    for i in range(NUMBER_OF_WORDS):
        if not crossword_words[i].visited:
            lst = []
            lst = get_list_of_component(i, crossword_words, lst)
            component_list.append(lst.copy())

    # make field 'visited' false, after DSF algorithm
    for i in range(NUMBER_OF_WORDS):
        crossword_words[i].visited = False

    return component_list


def put_positions(index, crossword_words: List[Word]):
    """
    This function determine position (horizontal or vertical) for each word
    This function uses DFS algorithm
    :param index: index of the word
    :param crossword_words: list of all words of the crosswords
    """
    if crossword_words[index].visited:
        return
    else:
        crossword_words[index].visited = True
    for nb in crossword_words[index].neighbours:
        # words which are crossed each other should have different positions
        if crossword_words[nb[1]].position == 1:
            if crossword_words[index].position != 0:
                crossword_words[index].position = 0
                break
        if crossword_words[nb[1]].position == 0:
            if crossword_words[index].position != 1:
                crossword_words[index].position = 1
                break

    if crossword_words[index].position == -1:
        # if no word has position from its 'neighbours', randomly choose position for this word
        crossword_words[index].position = random.randint(0, 1)

    # recursively call this function for all his 'neighbours'
    for nb in crossword_words[index].neighbours:
        put_positions(nb[1], crossword_words)


def put_word_on_grid(word: Word, grid: List[List[str]]):
    """
    This function puts word on a grid
    :param word: word of the crossword
    :param grid: given grid, where the word should be located
    :return: grid with placed word on it
    """
    global GRID_SIZE

    x = word.coordinates[0]
    y = word.coordinates[1]
    if word.position == 0:
        for c in word.name:
            grid[x][y] = c
            # safely increase coordinate of the cell
            x = (x + 1) % GRID_SIZE
    if word.position == 1:
        for c in word.name:
            grid[x][y] = c
            # safely increase coordinate of the cell
            y = (y + 1) % GRID_SIZE

    return grid


def put_connected_words_on_grid(index, crossword_words: List[Word], grid: List[List[str]]):
    """
    This function defines coordinates of a word
    and puts all words from one connectivity component on  a grid
    This function uses DFS algorithm
    :param index: index of a word from a connectivity component
    :param crossword_words: list of all words of this crossword
    :param grid: given frid
    :return: grid with the given word in it
    """
    global GRID_SIZE

    if crossword_words[index].visited:
        return grid
    else:
        crossword_words[index].visited = True

    for i, nb in enumerate(crossword_words[index].neighbours):
        # check if the word has already coordinates
        if crossword_words[nb[1]].coordinates[1] != -1:
            # coordinates of the word are uniquely determined form its 'neighbours'
            if crossword_words[index].position == 0:
                crossword_words[index].coordinates[0] = (crossword_words[nb[1]].coordinates[0]
                                                         - crossword_words[index].neighbours[i][0]) % GRID_SIZE
                crossword_words[index].coordinates[1] = (crossword_words[nb[1]].coordinates[1]
                                                         + crossword_words[index].neighbours[i][2]) % GRID_SIZE
            if crossword_words[index].position == 1:
                crossword_words[index].coordinates[0] = (crossword_words[nb[1]].coordinates[0]
                                                         + crossword_words[index].neighbours[i][2]) % GRID_SIZE
                crossword_words[index].coordinates[1] = (crossword_words[nb[1]].coordinates[1]
                                                         - crossword_words[index].neighbours[i][0]) % GRID_SIZE
            grid = put_word_on_grid(crossword_words[index], grid)
            break

    if crossword_words[index].coordinates[0] == -1:
        # if no words from 'neighbours' have coordinates, then put it in the center of a grid
        crossword_words[index].coordinates[0] = GRID_SIZE // 2
        crossword_words[index].coordinates[1] = GRID_SIZE // 2
        grid = put_word_on_grid(crossword_words[index], grid)

    # recursively call this function for all his 'neighbours'
    for nb in crossword_words[index].neighbours:
        grid = put_connected_words_on_grid(nb[1], crossword_words, grid)

    return grid


def create_list_of_grids(sample: List[Word]) -> List[List[List[str]]]:
    """
    This function creates a grid for all connectivity components of the crossword
    :param sample: list of all words of the crossword
    :return: List of all grids (one for every component connectivity)
    """
    global NUMBER_OF_WORDS
    global GRID_SIZE

    list_of_grids = []

    random_word_order = [i for i in range(NUMBER_OF_WORDS)]
    random.shuffle(random_word_order)

    for i in random_word_order:
        # go through every connectivity component
        if not sample[i].visited:
            picture = [['.' for _ in range(GRID_SIZE)] for _ in range(GRID_SIZE)]
            picture = put_connected_words_on_grid(i, sample, picture)
            list_of_grids.append(picture)

    # make field 'visited' false, after DSF algorithm
    for i in range(NUMBER_OF_WORDS):
        sample[i].visited = False

    return list_of_grids


def print_final_answer(crossword_words: List[Word], output_index: int):
    """
    This function prints final answer to output file
    :param crossword_words: list of all words of the crossword
    :param output_index: index of output file
    :return:
    """
    global GRID_SIZE
    global TARGET_GRID_SIZE
    global NUMBER_OF_WORDS

    x_min, y_min = GRID_SIZE - 1, GRID_SIZE - 1

    picture = [['.' for _ in range(TARGET_GRID_SIZE)] for _ in range(TARGET_GRID_SIZE)]

    # determine minimum x and y coordinates
    for wrd in crossword_words:
        x_min = min(x_min, wrd.coordinates[0])
        y_min = min(y_min, wrd.coordinates[1])

    # shift all coordinates so they can fit required size
    for i in range(NUMBER_OF_WORDS):
        crossword_words[i].coordinates = [crossword_words[i].coordinates[0] - x_min,
                                          crossword_words[i].coordinates[1] - y_min]

    # printing answer to output file
    with open(f'outputs/output{output_index}.txt', "w") as file1:
        for wrd in crossword_words:
            # change position of words, because in the code 'horizontal' words lies along OX and 'vertical' along OY
            if wrd.position == 0:
                pos = 1
            else:
                pos = 0
            file1.write(f'{wrd.coordinates[0]} {wrd.coordinates[1]} {pos}\n')

    # uncomment if you need graphical representation of a crossword in console

    for wrd in crossword_words:
        x = wrd.coordinates[0]
        y = wrd.coordinates[1]
        for c in wrd.name:
            picture[x][y] = c
            if wrd.position == 0:
                x += 1
            else:
                y += 1

    print(f'Solution of {output_index} crossword')

    for x in range(TARGET_GRID_SIZE):
        for y in range(TARGET_GRID_SIZE):
            print(picture[x][y], end=" ")
        print()


def delete_useless_letter(crossword_dictionary: Dict[str, List[int]]):
    """
    This function deletes letters, which are only in ine word, or are not in any word
    :param crossword_dictionary: crossword dictionary
    :return: without useless letters
    """

    useless_keys = [key for key, value in crossword_dictionary.items() if
                    (not value) or len(value) == 1 or identical_elements(value)]
    for key in useless_keys:
        del crossword_dictionary[key]

    return crossword_dictionary


def link_crossword_words(sample_indexes_original: List[int], sample_dictionary: Dict[str, List[int]],
                         crossword_words_original: List[Word]) -> Tuple[List[Word], Dict[str, List[int]]]:
    """
    This function creates random crossword, based on one requirement: every word should be connected at least to another
    one, if it is possible, but it can't connect words, which already connected with another one
    :param sample_indexes_original: list of words' indexes, which have already a 'neighbour'
    :param sample_dictionary: dictionary of current crossword
    :param crossword_words_original: list of all crossword words
    :return: list of all words and dictionary of new crossword
    """
    global words
    global NUMBER_OF_WORDS

    # coping the info to avoid changes in input data
    crossword_dictionary = copy.deepcopy(sample_dictionary)
    tmp_dict = copy.deepcopy(crossword_dictionary)
    crossword_words = copy.deepcopy(crossword_words_original)
    sample_indexes = copy.deepcopy(sample_indexes_original)

    # delete useless letters from the dictionary
    tmp_dict = delete_useless_letter(tmp_dict)

    while True:
        # check number of words which have at least one 'neighbour' and check if we can connect any other words
        if len(sample_indexes) == NUMBER_OF_WORDS or (not tmp_dict):
            break
        # randomly choose letter, by which we will connect two words
        random_key = random.choice(list(tmp_dict.keys()))
        # copy indexes of all elements, which have this letter 'free'
        random_indexes = tmp_dict[random_key]
        # this list will have indexes of words which has already at least one 'neighbour'
        exist_array = []
        # this list will have indexes of words which don't have any 'neighbour'
        not_exist_array = []
        for idx in random_indexes:
            if idx in sample_indexes:
                exist_array.append(idx)
            else:
                not_exist_array.append(idx)

        if len(exist_array) == len(random_indexes):
            # all words having the chosen letter have already at least one 'neighbour',
            # therefore we can delete this letter
            del tmp_dict[random_key]
            continue

        if len(exist_array) > 0:
            # randomly choose one word from each array
            idx1 = random.choice(exist_array)
            idx2 = random.choice(not_exist_array)
            sample_indexes.append(idx2)
        else:
            # randomly chose two different words from this letter
            indexes = random.sample(list(set(random_indexes)), 2, counts=None)
            idx1 = indexes[0]
            idx2 = indexes[1]
            # add them to the list which contains indexes of all words, which have at least one 'neighbour'
            sample_indexes.append(idx1)
            sample_indexes.append(idx2)

        # remove the indexes from the key list
        crossword_dictionary[random_key].remove(idx1)
        crossword_dictionary[random_key].remove(idx2)
        tmp_dict[random_key].remove(idx1)
        tmp_dict[random_key].remove(idx2)

        # if the word has several 'random' letters, randomly choose one of them
        char_index_array = []
        for idx, c in enumerate(words[idx1]):
            if c == random_key:
                char_index_array.append(idx)
        for c in crossword_words[idx1].neighbours:
            # check if the letter is already busy
            if c[0] in char_index_array:
                char_index_array.remove(c[0])
        idx_of_c1 = random.choice(char_index_array)

        # if the word has several 'random' letters, randomly choose one of them
        char_index_array = []
        for idx, c in enumerate(words[idx2]):
            if c == random_key:
                char_index_array.append(idx)
        idx_of_c2 = random.choice(char_index_array)

        # add chosen words as neighbours to each word
        crossword_words[idx1].neighbours.append([idx_of_c1, idx2, idx_of_c2])
        crossword_words[idx2].neighbours.append([idx_of_c2, idx1, idx_of_c1])

        # clean up the dictionary
        tmp_dict = delete_useless_letter(tmp_dict)

        # check number of connected words
        if len(sample_indexes) == NUMBER_OF_WORDS:
            break

    # put positions for each word
    for i in range(NUMBER_OF_WORDS):
        if (not crossword_words[i].visited) and (crossword_words[i].position != -1):
            put_positions(i, crossword_words)

    for i in range(NUMBER_OF_WORDS):
        crossword_words[i].visited = False

    # put positions for each word
    for i in range(NUMBER_OF_WORDS):
        if not crossword_words[i].visited:
            put_positions(i, crossword_words)

    for i in range(NUMBER_OF_WORDS):
        crossword_words[i].visited = False

    return crossword_words, crossword_dictionary


def mutation(crossword_words: List[Word], crossword_dict: Dict[str, List[int]]):
    """
    This function implements mutation of the crossword - detach one word from everyone and reattach it
    :param crossword_words: list of al crossword words
    :param crossword_dict: dictionary of the crossword
    :return: mutated dictionary
    """
    global NUMBER_OF_WORDS

    # coping the info to avoid changes in input data
    mutation_dict = copy.deepcopy(crossword_dict)
    mutation_words = copy.deepcopy(crossword_words)

    # randomly choose index of a word which we want to reattach
    index1 = random.randint(0, NUMBER_OF_WORDS - 1)

    # detach chosen word from all his neighbours
    for nb in mutation_words[index1].neighbours:
        letter = mutation_words[index1].name[nb[0]]
        mutation_dict[letter].append(index1)
        mutation_dict[letter].append(nb[1])
        mutation_words[nb[1]].neighbours = (
            list(filter(lambda x: x[1] != index1, mutation_words[nb[1]].neighbours)))
    mutation_words[index1].neighbours = []

    # get list of all components
    component_array = get_list_of_all_components(mutation_words)
    component_array = sorted(component_array, key=lambda x: len(x), reverse=True)

    # determine indexes of all words which have at least one neighbour
    sample_indexes = []
    for cmp in component_array:
        if len(cmp) > 1:
            for i in cmp:
                sample_indexes.append(i)

    # reset all unconnected words
    for i in range(NUMBER_OF_WORDS):
        if not (i in component_array[0]):
            mutation_words[i].position = -1
        mutation_words[i].coordinates = [-1, -1]
        mutation_words[i].visited = False

    # try to connect unconnected words
    return link_crossword_words(sample_indexes, mutation_dict, mutation_words)


def crossover(crossover_words_1: List[Word], crossover_words_2: List[Word]):
    """
    This function implements crossover
    :param crossover_words_1: list of  the 1-st crossword words
    :param crossover_words_2: list of the 2-nd crossword words
    :return: list of new crossword words and its dictionary
    """
    global alphabet_dict
    global NUMBER_OF_WORDS
    global words

    # coping the info to avoid changes in input data
    child_dict = copy.deepcopy(alphabet_dict)
    mother_words = copy.deepcopy(crossover_words_1)
    father_words = copy.deepcopy(crossover_words_2)

    # number of words, which are taken from mother crossword
    mother_number = NUMBER_OF_WORDS // 2

    # get list of all components from mother crossword
    mother_components = get_list_of_all_components(mother_words)
    mother_components = sorted(mother_components, key=lambda x: len(x), reverse=True)

    # randomly choose component of length 'mother_number' from mother crossword
    if len(mother_components[0]) >= mother_number:
        mother_components = list(filter(lambda x: len(x) >= mother_number, mother_components))
        mother_indexes = mother_components[random.randint(0, len(mother_components) - 1)]
        starting_index = random.choice(mother_indexes)
        mother_indexes = []
        mother_indexes = get_list_of_component(starting_index, mother_words, mother_indexes)
        # write indexes of all mother's words
        mother_indexes = mother_indexes[:mother_number]
    else:
        # write indexes of all mother's words
        mother_indexes = mother_components[0]

    # write indexes of all father's words
    father_indexes = []
    for i in range(NUMBER_OF_WORDS):
        if not (i in mother_indexes):
            father_indexes.append(i)

    # create list of words for child crossword
    child_words = [Word(w) for w in words]

    # connect words from 'mother_indexes' in child crossword as they were connected in mother crossword
    tmp_indexes = []
    for i in mother_indexes:
        tmp_indexes.append(i)
        child_words[i].position = mother_words[i].position
        for nb in mother_words[i].neighbours:
            if not (nb[1] in tmp_indexes) and nb[1] in mother_indexes:
                letter = mother_words[i].name[nb[0]]
                child_dict[letter].remove(i)
                child_dict[letter].remove(nb[1])
                child_words[i].neighbours.append(nb)
                child_words[nb[1]].neighbours.append([nb[2], i, nb[0]])

    # connect words from 'father_indexes' in child crossword as they were connected in father crossword
    tmp_indexes = []
    for i in father_indexes:
        tmp_indexes.append(i)
        for nb in father_words[i].neighbours:
            if not (nb[1] in tmp_indexes) and nb[1] in father_indexes:
                letter = father_words[i].name[nb[0]]
                child_dict[letter].remove(i)
                child_dict[letter].remove(nb[1])
                child_words[i].neighbours.append(nb)
                child_words[nb[1]].neighbours.append([nb[2], i, nb[0]])

    # get list of all components in child crossword
    crossover_components = get_list_of_all_components(child_words)

    # remove component containing mother component
    for cmp in crossover_components:
        if cmp[0] in mother_indexes:
            crossover_components.remove(cmp)
            break

    random.shuffle(crossover_components)

    # trying to link each component to mother component, as they were connected in father crossword
    while crossover_components:
        component = crossover_components[-1]
        random.shuffle(component)

        ok = True
        for idx in component:
            if not ok:
                break
            for nb in father_words[idx].neighbours:
                if not ok:
                    break
                if not (nb[1] in father_indexes):
                    flag = True
                    for nbnb in child_words[nb[1]].neighbours:
                        if nbnb[0] == nb[2]:
                            flag = False
                            break
                    if flag:
                        letter = child_words[idx].name[nb[0]]
                        child_dict[letter].remove(idx)
                        child_dict[letter].remove(nb[1])
                        child_words[idx].neighbours.append(nb)
                        child_words[nb[1]].neighbours.append([nb[2], idx, nb[0]])
                        ok = False
                        break

        crossover_components.pop()

    # get list of all components of child crossword
    crossover_components = get_list_of_all_components(child_words)
    # indexes of all words which don't have a 'neighbour'
    sample_indexes = []

    for cmp in crossover_components:
        if len(cmp) > 1:
            for i in cmp:
                sample_indexes.append(i)

    # try to connect unconnected words in child crossword
    return link_crossword_words(sample_indexes, child_dict, child_words)


def indicate_hard_case(x1, y1, x2, y2, word_list: List[Word], grid: List[List[str]]):
    """
    This function returns True if words have the same 'position' and have one letter near each other
    :param x1:
    :param y1:
    :param x2:
    :param y2:
    :param word_list: list of all crossword word
    :param grid: grid of the crossword
    :return: True if the words are staying in incorrect positions, otherwise - False
    """
    global GRID_SIZE

    if x1 != x2:
        # both words are vertical
        y = y1
        x = x1
        counter = 0
        while grid[x][y] != '.':
            x = (x - 1) % GRID_SIZE
            counter += 1
            if counter >= GRID_SIZE:
                return True
        wrd_string = ""
        x = (x + 1) % GRID_SIZE
        x_start, y_start = x, y
        # find a word which is created by theirs 'connection'
        while grid[x][y] != '.':
            wrd_string += grid[x][y]
            x = (x + 1) % GRID_SIZE
        # check if this word exists in the crossword
        for wrd in word_list:
            if wrd.name == wrd_string and wrd.coordinates == [x_start, y_start] and wrd.position == 0:
                return False
    if y1 != y2:
        # both words are horizontal
        x = x1
        y = y1
        counter = 0
        while grid[x][y] != '.':
            y = (y - 1) % GRID_SIZE
            counter += 1
            if counter >= GRID_SIZE:
                return True
        wrd_string = ""
        y = (y + 1) % GRID_SIZE
        x_start, y_start = x, y
        # find a word which is created by theirs 'connection'
        while grid[x][y] != '.':
            wrd_string += grid[x][y]
            y = (y + 1) % GRID_SIZE
        # check if this word exists in the crossword
        for wrd in word_list:
            if wrd.name == wrd_string and wrd.coordinates == [x_start, y_start] and wrd.position == 1:
                return False
    return True


def fitness_function(crossword_words: List[Word], grid_list: List[List[List[str]]]):
    """
    This function calculates fitness function of the crossword, based on its grid
    :param crossword_words: list of all crossword words
    :param grid_list: list of all grids (components) of the crossword
    :return: numbers of errors in the crossword
    """
    global NUMBER_OF_WORDS
    global GRID_SIZE
    global TARGET_GRID_SIZE

    # number of errors in a crossword
    errors = 0

    # penalty coefficient for each error
    connectivity_component_penalty = -10
    wrong_neighbour_penalty = -1
    wrong_word_cross_penalty = -1
    wrong_size_penalty = -1

    # minimum and maximum coordinates of all words
    x_max = 0
    y_max = 0
    x_min = GRID_SIZE - 1
    y_min = GRID_SIZE - 1

    # get list of all components
    components = get_list_of_all_components(crossword_words)

    # check number of components
    if len(components) > 1:
        errors += (len(components) - 1) * connectivity_component_penalty

    # for each component check number of errors
    for component_index, cmp in enumerate(components):
        # get grid of the component
        component_grid = grid_list[component_index]
        # list containing starting and end letter coordinates of the word
        start_end_letter_indexes = []
        for wrd_index in cmp:
            # get coordinates of the first letter of the word
            x = crossword_words[wrd_index].coordinates[0]
            y = crossword_words[wrd_index].coordinates[1]
            start_end_letter_indexes.append([x, y])
            x_min, x_max = (min(x_min, x), max(x_max, x))
            y_min, y_max = (min(y_min, y), max(y_max, y))
            # indicates the number of consecutive characters adjacent to the word on one side
            neighbour_number_1 = 0
            # indicates the number of consecutive characters adjacent to the word on another side
            neighbour_number_2 = 0
            if crossword_words[wrd_index].position == 0:
                # word has horizontal position

                # check if the word has another letter above it
                if component_grid[(x - 1) % GRID_SIZE][y] != '.':
                    errors += wrong_neighbour_penalty

                x = (x - 1) % GRID_SIZE
                # check correctness of the word on the grid
                for c in crossword_words[wrd_index].name:
                    x = (x + 1) % GRID_SIZE
                    if component_grid[x][y] != c:
                        errors += wrong_word_cross_penalty

                    # check adjacent character on one side
                    if component_grid[x][(y + 1) % GRID_SIZE] != '.':
                        neighbour_number_1 += 1
                        if [x, (y + 1) % GRID_SIZE] in start_end_letter_indexes:
                            if indicate_hard_case(x, y, x, (y + 1) % GRID_SIZE, crossword_words,
                                                  component_grid):
                                errors += wrong_neighbour_penalty
                                neighbour_number_1 = 0
                    else:
                        neighbour_number_1 = 0

                    # check adjacent character on another side
                    if component_grid[x][(y - 1) % GRID_SIZE] != '.':
                        neighbour_number_2 += 1
                        if [x, (y - 1) % GRID_SIZE] in start_end_letter_indexes:
                            if indicate_hard_case(x, y, x, (y - 1) % GRID_SIZE, crossword_words,
                                                  component_grid):
                                errors += wrong_neighbour_penalty
                                neighbour_number_2 = 0
                    else:
                        neighbour_number_2 = 0

                    # penalty if the word has consecutive adjacent characters
                    if neighbour_number_1 > 1 or neighbour_number_2 > 1:
                        errors += wrong_neighbour_penalty

                    x_min, x_max = (min(x_min, x), max(x_max, x))

                start_end_letter_indexes.append([x, y])

                # check if the word has another letter below it
                if component_grid[(x + 1) % GRID_SIZE][y] != '.':
                    errors += wrong_neighbour_penalty
            else:
                # word has vertical position

                # check if the word has another letter to its left
                if component_grid[x][(y - 1) % GRID_SIZE] != '.':
                    errors += wrong_neighbour_penalty

                y = (y - 1) % GRID_SIZE
                # check correctness of the word on the grid
                for c in crossword_words[wrd_index].name:
                    y = (y + 1) % GRID_SIZE
                    if component_grid[x][y] != c:
                        errors += wrong_word_cross_penalty

                    # check adjacent character on one side
                    if component_grid[(x + 1) % GRID_SIZE][y] != '.':
                        neighbour_number_1 += 1
                        if [(x + 1) % GRID_SIZE, y] in start_end_letter_indexes:
                            if indicate_hard_case(x, y, (x + 1) % GRID_SIZE, y, crossword_words,
                                                  component_grid):
                                errors += wrong_neighbour_penalty
                                neighbour_number_1 = 0
                    else:
                        neighbour_number_1 = 0

                    # check adjacent character on another side
                    if component_grid[(x - 1) % GRID_SIZE][y] != '.':
                        neighbour_number_2 += 1
                        if [(x - 1) % GRID_SIZE, y] in start_end_letter_indexes:
                            if indicate_hard_case(x, y, (x - 1) % GRID_SIZE, y, crossword_words,
                                                  component_grid):
                                errors += wrong_neighbour_penalty
                                neighbour_number_2 = 0
                    else:
                        neighbour_number_2 = 0

                    if neighbour_number_1 > 1 or neighbour_number_2 > 1:
                        errors += wrong_neighbour_penalty

                    y_min, y_max = (min(y_min, y), max(y_max, y))
                start_end_letter_indexes.append([x, y])

                # check if the word has another letter to its right
                if component_grid[x][(y + 1) % GRID_SIZE] != '.':
                    errors += wrong_neighbour_penalty

    # check size of the crossword
    if abs(x_max - x_min) + 1 > TARGET_GRID_SIZE:
        errors += (abs(x_max - x_min) - TARGET_GRID_SIZE + 1) * wrong_size_penalty
    if abs(y_max - y_min) + 1 > TARGET_GRID_SIZE:
        errors += (abs(y_max - y_min) - TARGET_GRID_SIZE + 1) * wrong_size_penalty

    return errors


def create_initial_population() -> List[Sample]:
    """
    This function creates initial populations of crosswords
    :return: list of crosswords
    """
    global INITIAL_POPULATION_NUMBER
    global POPULATION_NUMBER
    global words

    # list of all generated crosswords
    samples = []
    # list of all crosswords words
    tmp_list_of_words = [Word(w) for w in words]
    # list of all connected words - initially empty
    empty_list = []
    # creating population
    for i in range(INITIAL_POPULATION_NUMBER):
        wrd_array, smp_dict = link_crossword_words(empty_list, alphabet_dict, tmp_list_of_words)
        grd_arr = create_list_of_grids(wrd_array)
        ft_fnc = fitness_function(wrd_array, grd_arr)
        samples.append(Sample(wrd_array, smp_dict, grd_arr, ft_fnc))

    # sort all crosswords by their fitness function
    samples = sorted(samples, key=lambda x: x.fitness_function_value, reverse=True)
    # leave only the strongest POPULATION_NUMBER crosswords
    samples = list(samples[:POPULATION_NUMBER])

    return samples


def make_new_generation(samples: List[Sample]) -> List[Sample]:
    """
    This function creates new generation of the population
    :param samples: list of crosswords (population)
    :return: list of crosswords (population)
    """
    global MUTATION_RATE
    global OFFSPRING_RATE
    global POPULATION_NUMBER

    # calculate number of mutations for the population
    mutation_number = int(POPULATION_NUMBER * MUTATION_RATE)
    # randomly choose mutation_number crosswords from population
    mutation_indexes = random.sample(range(len(samples)), mutation_number, counts=None)
    # mutate chosen crosswords
    for i in mutation_indexes:
        wrd_array, smp_dict = mutation(samples[i].word_list, samples[i].dictionary)
        pct_arr = create_list_of_grids(wrd_array)
        ft_fnc = fitness_function(wrd_array, pct_arr)
        samples.append(Sample(wrd_array, smp_dict, pct_arr, ft_fnc))

    # calculate number of new 'children'
    crossover_number = int(POPULATION_NUMBER * OFFSPRING_RATE)
    # randomly choose crossover_number crosswords from population
    crossover_indexes = [random.sample(range(len(samples)), 2) for _ in range(crossover_number)]
    # crossover chosen crosswords
    for sample_pair in crossover_indexes:
        wrd_array, smp_dict = crossover(samples[sample_pair[0]].word_list, samples[sample_pair[1]].word_list)
        pct_arr = create_list_of_grids(wrd_array)
        ft_fnc = fitness_function(wrd_array, pct_arr)
        samples.append(Sample(wrd_array, smp_dict, pct_arr, ft_fnc))

    random.shuffle(samples)
    # sort all crosswords by their fitness function
    samples = sorted(samples, key=lambda x: x.fitness_function_value, reverse=True)
    # leave only the strongest POPULATION_NUMBER crosswords
    samples = list(samples[:POPULATION_NUMBER])

    return samples


def make_crossword(crossword_index: int):
    """
    This function makes valid crossword
    :param crossword_index: number of crossword (input file)
    """
    global NUMBER_OF_WORDS
    global alphabet_dict
    global MAXIMUM_NUMBER_OF_SAME_FITNESS_FUNCTION
    global words
    global POPULATION_NUMBER
    global INITIAL_POPULATION_NUMBER

    # maximum number of attempts make a crossword
    maximum_attempts_number = 3

    # list containing all words
    words = []
    # read words from input file
    with open(f'inputs/input{crossword_index}.txt', 'r') as file:
        # Read lines and remove leading/trailing whitespaces
        words = [line.strip() for line in file]

    # change number of wards variable
    NUMBER_OF_WORDS = len(words)

    # for each letter (key) add index of the word which contains it
    alphabet_dict = {chr(i): [] for i in range(ord('a'), ord('z') + 1)}
    for word_index, word in enumerate(words):
        for letter in word:
            alphabet_dict[letter].append(word_index)

    # delete useless letters from dictionary
    alphabet_dict = delete_useless_letter(alphabet_dict)

    # create initial population
    samples = create_initial_population()

    # number of generation
    generation = 0
    # current number of attempts make a crossword
    attempt_counter = 0
    # max previous fitness function of previous generation
    prev_max_fitness_function = samples[0].fitness_function_value
    # show how many generations was the same max fitness function
    same_fitness_function_counter = 0
    while samples[0].fitness_function_value != 0:
        # Uncomment if you need information about current generation
        print(
            f'generation = {generation}, max fitness function = {samples[0].fitness_function_value}, # {crossword_index}')
        samples = make_new_generation(samples)
        generation += 1
        # check if max fitness function was changed
        if samples[0].fitness_function_value == prev_max_fitness_function:
            same_fitness_function_counter += 1
        else:
            same_fitness_function_counter = 0
            prev_max_fitness_function = samples[0].fitness_function_value

        if same_fitness_function_counter > MAXIMUM_NUMBER_OF_SAME_FITNESS_FUNCTION:
            # it means that we are in local maximum, and we need to renew thr population
            POPULATION_NUMBER = 100
            INITIAL_POPULATION_NUMBER = 1000
            samples = create_initial_population()
            generation = 0
            attempt_counter += 1
            prev_max_fitness_function = samples[0].fitness_function_value
            same_fitness_function_counter = 0
            if attempt_counter > maximum_attempts_number:
                # we didn't manage to create a crossword
                with open(f'outputs/output{crossword_index}.txt', "w") as file1:
                    file1.write("\n")
                return

    # printing answer to the file
    print_final_answer(samples[0].word_list, crossword_index)


# number of words in crossword
NUMBER_OF_WORDS = 0
# number of crosswords in population
POPULATION_NUMBER = 100
# number of crosswords in initial population
INITIAL_POPULATION_NUMBER = 1000
# percentage of 'POPULATION_NUMBER' of offsprings after mutation,
MUTATION_RATE = 1
# percentage of 'POPULATION_NUMBER' of offsprings after crossover,
OFFSPRING_RATE = 1
# maximum number of generation till we reach local maximum
MAXIMUM_NUMBER_OF_SAME_FITNESS_FUNCTION = 200
# size of a grid where crosswords are written
GRID_SIZE = 50
# size of a grid where we need to put all words
TARGET_GRID_SIZE = 20

# dictionary, where key - all lower characters of alphabet
alphabet_dict = {chr(i): [] for i in range(ord('a'), ord('z') + 1)}
# initial list of words
words = []

# generate 'NUMBER_OF_CROSSWORDS' crosswords
for i in range(1, NUMBER_OF_CROSSWORDS + 1):
    try:
        make_crossword(i)
    except Exception as e:
        continue
