import time
import re


def contatenate_path(paths, tagLen, pattern):
    '''path = [[[1,4,5,6], 'ADR'], [], []]'''
    res = dict()
    for idx1 in range(len(paths)):
        if len(paths[idx1][0]) > tagLen:
            flag = True  # add the result at the end
            for idx2 in range(len(paths)):
                if len(paths[idx2][0]) > tagLen:
                    if paths[idx1][0] == paths[idx2][0]:  # same tags
                        continue
                    path1 = paths[idx1]  # shorter length
                    path2 = paths[idx2]
                    if len(path2[0]) < len(path1[0]):
                        path1, path2 = path2, path1
                    for n in range(len(path1[0]) + len(path2[0]) - 1):
                        x1 = max(len(path1[0]) - 1 - n, 0)
                        y1 = min(len(path1[0]) + len(path2[0]) - n - 1, len(path1[0]))
                        x2 = max(n + 1 - len(path1[0]), 0)
                        y2 = min(n + 1, len(path2[0]))
                        if path1[0][x1:y1] == path2[0][x2:y2]:
                            result1 = re.findall(pattern, path1[1])
                            result2 = re.findall(pattern, path2[1])
                            if n < len(path1[0]) - 1:
                                new_tag = path1[0] + path2[0][n+1:]
                                new_str = result1 + result2[n:]
                            elif n > len(path2[0]) - 1:
                                new_tag = path2[0] + path1[0][len(path1[0]) + len(path2[0]) - 1 - n:]
                                new_str = result2 + result1[len(path1[0]) + len(path2[0]) - 2 - n:]
                            else:
                                new_tag = path2[0]
                                new_str = result2
                            res.setdefault('-'.join([str(ele) for ele in new_tag]), set())
                            res['-'.join([str(ele) for ele in new_tag])].add(''.join(new_str))
                            flag = False
                            break
            if flag:
                res.setdefault('-'.join([str(ele) for ele in paths[idx1][0]]), set())
                res['-'.join([str(ele) for ele in paths[idx1][0]])].add(paths[idx1][1])

    new_res = []
    for i, j in res.items():
        ls = [int(ele) for ele in i.split('-')]
        for jj in j:
            new_res.append((ls, jj))
    return new_res


def find_all_paths(graph, start):
    paths = []
    endpoints = set()  # Set to store encountered endpoints
    visited = set()  # Set to keep track of visited nodes

    def dfs(current_node, path, current_weight):
        visited.add(current_node)
        path.append(current_node)
        if graph.get(current_node) is None:
            graph[current_node] = []
        if len(graph[current_node]) == 0:
            endpoints.add(current_node)  # Record endpoint

        for neighbor, weight in graph[current_node]:
            if neighbor not in visited:
                dfs(neighbor, path, current_weight + weight)

        visited.remove(current_node)
        path.pop()

    dfs(start, [], '')

    for endpoint in endpoints:
        endpoint_paths = find_paths_to_endpoint(graph, start, endpoint)
        paths.extend(endpoint_paths)

    return paths


def find_paths_to_endpoint(graph, start, endpoint):
    paths = []
    visited = set()  # Set to keep track of visited nodes

    def dfs(current_node, path, current_weight):
        visited.add(current_node)
        path.append(current_node)

        if current_node == endpoint:
            paths.append((path[:], current_weight))  # Append a copy of the current path and its weight

        for neighbor, weight in graph[current_node]:
            if neighbor not in visited:
                dfs(neighbor, path, current_weight + weight)

        visited.remove(current_node)
        path.pop()

    dfs(start, [], '')
    return paths


def graph_path_filter(paths, specIntensity, topK=10, tagLen=4):
    """return top k different tag results, sorted by tag length and total intensity,
     one series of tags could represent multiple sequences"""
    filter_paths = list()  # first check if any peak is extremely low
    extreme_factor = 10.0
    pattern = r'\([A-Za-z]+\)|[A-Za-z]'
    for path in paths:
        if len(path[0]) > tagLen and path[1].count('(') <= 1:  # at most one discontinued AAs in the tag:
            break_points = list()
            for idx in range(1, len(path[0]) - 1):
                if extreme_factor * specIntensity[path[0][idx]] < \
                        (specIntensity[path[0][idx - 1]] + specIntensity[path[0][idx + 1]]) / 2:
                    break_points.append(idx)
            if len(break_points) == 0:
                filter_paths.append(path)
            if len(break_points) == 1:
                result = re.findall(pattern, path[1])
                if ''.join(result) != path[1]:
                    raise Exception("Sorry, invalid modification format exists!")
                if break_points[0] > tagLen:
                    filter_paths.append((path[0][:break_points[0]], ''.join(result[:break_points[0] - 1])))
                if len(path[0]) - break_points[0] > tagLen + 1:
                    pass
    common_factor_paths = dict()  # check similarity and use common factor
    common_threshold = 0.7
    common_base = 500
    filter_paths.sort(key=lambda s: (-len(s[0]), -sum([specIntensity[ele1] for ele1 in s[0]])))
    for path in filter_paths[:common_base]:
        flag = False  # similar ones
        for reference in filter_paths[:common_base]:
            if path[0] == reference[0] or len(set(path[0]).intersection(set(reference[0]))) / len(path[0]) == 1:  # same path
                continue
            elif len(set(path[0]).intersection(set(reference[0]))) / len(path[0]) > common_threshold:
                flag = True
                result = re.findall(pattern, path[1])
                diff_list = sorted(set(path[0]) - set(reference[0]))
                # print(path[0])
                # print(reference[0])
                # print(diff_list)
                # print('next')
                diff_idx = [path[0].index(num) for num in diff_list]
                if diff_idx[0] - 1 >= tagLen:
                    key = '-'.join([str(num) for num in path[0][0: diff_idx[0]]])
                    value = ''.join(result[0: diff_idx[0] - 1])
                    common_factor_paths.setdefault(key, set())
                    common_factor_paths[key].add(value)
                if len(path[0]) - diff_idx[-1] - 2 >= tagLen:
                    key = '-'.join([str(num) for num in path[0][diff_idx[-1] + 1:]])
                    value = ''.join(result[diff_idx[-1] + 1:])
                    common_factor_paths.setdefault(key, set())
                    common_factor_paths[key].add(value)
                if len(diff_idx) >= 2:
                    for idx in range(len(diff_idx) - 1):
                        if diff_idx[idx + 1] - diff_idx[idx] - 2 >= tagLen:
                            key = '-'.join([str(num) for num in path[0][diff_idx[idx] + 1: diff_idx[idx + 1]]])
                            value = ''.join(result[diff_idx[idx] + 1: diff_idx[idx + 1] - 1])
                            common_factor_paths.setdefault(key, set())
                            common_factor_paths[key].add(value)
        if not flag:
            key = '-'.join([str(num) for num in path[0]])
            value = path[1]
            common_factor_paths.setdefault(key, set())
            common_factor_paths[key].add(value)
    filter_paths = list()
    for k, v in common_factor_paths.items():
        for ele in v:
            path = [int(i) for i in k.split('-')]
            filter_paths.append((path, ele))
    '''concatenate two if possible'''
    # filter_paths = contatenate_path(filter_paths, tagLen, pattern)
    # for i in sorted(filter_paths):
    #     print('after', i)

    filter_paths.sort(key=lambda s: (-len(s[0]), -sum([specIntensity[ele] for ele in s[0]])))
    # for i in filter_paths:
    #     print(i[0], i[1], [round(specIntensity[ii], 2) for ii in i[0]])
    unique = set()  # top k tags
    res = list()
    for path in filter_paths:
        mismatch = 1 if len(path[0]) > tagLen + 1 else 0
        if len(path[0]) < tagLen + 1:  # peak number
            break
        unique.add('-'.join([str(x) for x in path[0]]))
        if len(unique) <= topK:

            simbol = path[1].replace('(', '').replace(')', '')
            if 0 in path[0]:  # should be b ion
                res.append([path[0], path[1], simbol, mismatch])
            elif 1 in path[0]:  # should be y ion
                res.append([path[0], path[1], simbol[::-1], mismatch])
            else:
                res.append([path[0], path[1], simbol, mismatch])
                res.append([path[0], path[1], simbol[::-1], mismatch])
        # if len(unique) <= topK and path[1].count('(') <= 2:  # at most two discontinued AAs in the tag
        #     pos = [p for p in range(len(path[1])) if path[1][p] == '(']
        #
        #     mismatch = 1
        #     valid_len = len(path[1]) - 2 * len(pos)
        #     if path[1].count('(') == 2:
        #         mismatch = 2
        #     elif path[1].count('(') == 1 and valid_len > tagLen + 1:
        #         mismatch = 2
        #     elif valid_len <= tagLen:
        #         mismatch = 0
        #
        #     simbol = list(path[1])
        #     for p in pos:
        #         if simbol[p + 1] > simbol[p + 2]:
        #             simbol[p + 2] = simbol[p + 1]
        #         else:
        #             simbol[p + 1] = simbol[p + 2]
        #     simbol = ''.join([sim for sim in simbol if sim not in ['(', ')']])
        #     res.append([path[0], path[1], simbol, mismatch])
        #     res.append([path[0], path[1], simbol[::-1], mismatch])
    return res





if __name__ == "__main__":

    # Example graph
    graph1 = {
        1: [(2, 'C'), (3, 'A')],
        2: [(4, 'I'), (4, 'L')],
        3: [(4, 'S'), (5, 'T'), (8, 'D')],
        4: [(5, 'V'), (6, 'N'), (6, '#')],
        5: [(7, 'M')],
        9: [(10, 'E')]

    }
    # graph = {
    #     'A': [('B', 2), ('B',2.5), ('B', 1)],
    #     'B': []
    # }

    start = 111
    t1 = time.perf_counter()
    paths1 = find_all_paths(graph1, start)
    t2 = time.perf_counter()
    print(t2 - t1)
    print(paths1)
    print(graph1)

    '''####################################'''

