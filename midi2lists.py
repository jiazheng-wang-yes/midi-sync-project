import sys
import mido
from mido import MidiFile, MetaMessage, bpm2tempo
import numpy as np

def get_notes_with_start_time(midi_file_path):
    midi_file = MidiFile(midi_file_path)
    notes = []

    for track in midi_file.tracks:
        current_time = 0

        for msg in track:
            current_time += msg.time

            if msg.type == 'note_on' and msg.velocity > 0:
                notes.append((msg.note, current_time))

    return sorted(notes, key = lambda x : x[1])
def filter_tuples(tuples_list):
    filtered_dict = {}
    
    for t in tuples_list:
        if t[1] not in filtered_dict:
            filtered_dict[t[1]] = t[0]
        else:
            filtered_dict[t[1]] = max(filtered_dict[t[1]], t[0])
    tuples_list = [(value, key) for key, value in filtered_dict.items()]

    
    return tuples_list
def find_local_maximums(tuple_list):
    tuple_list.insert(0, (float('-inf'), 0))
    tuple_list.append((float('-inf'), float('inf')))
    local_maxima = []

    for i in range(1, len(tuple_list) - 1):
        if tuple_list[i][0] > tuple_list[i - 1][0] and tuple_list[i][0] >= tuple_list[i + 1][0]:
            local_maxima.append(tuple_list[i])

    return local_maxima
def needleman_wunsch(seq1, seq2, gap_penalty=-1, match_score=1, mismatch_penalty=-1):
    m, n = len(seq1), len(seq2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    traceback = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(1, m + 1):
        dp[i][0] = gap_penalty * i
        traceback[i][0] = "up"

    for j in range(1, n + 1):
        dp[0][j] = gap_penalty * j
        traceback[0][j] = "left"

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = dp[i - 1][j] + gap_penalty
            insert = dp[i][j - 1] + gap_penalty

            dp[i][j] = max(match, delete, insert)

            if dp[i][j] == match:
                traceback[i][j] = "down"
            elif dp[i][j] == delete:
                traceback[i][j] = "up"
            else:
                traceback[i][j] = "left"

    aligned_seq1, aligned_seq2 = [], []
    i, j = m, n
    while i > 0 or j > 0:
        if traceback[i][j] == "down":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif traceback[i][j] == "up":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(None)
            i -= 1
        else:
            aligned_seq1.append(None)
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    return aligned_seq1[::-1], aligned_seq2[::-1]


def match_pitches_and_return_times(list1, list2):
    pitches1, times1 = zip(*list1)
    pitches2, times2 = zip(*list2)

    aligned_pitches1, aligned_pitches2 = needleman_wunsch(list(pitches1), list(pitches2))

    aligned_times1, aligned_times2 = [], []
    i1, i2 = 0, 0
    for p1, p2 in zip(aligned_pitches1, aligned_pitches2):

        if p1 is not None:
            aligned_times1.append(times1[i1])
            i1 += 1
        else:
            aligned_times1.append(None)

        if p2 is not None:
            aligned_times2.append(times2[i2])
            i2 += 1
        else:
            aligned_times2.append(None)
    output1, output2 = [], []
    for a, b in zip(aligned_times1, aligned_times2):
        if (a is not None and b is not None):
            output1.append(a)
            output2.append(b)
    return output1, output2


def check_point_mapping(a, b):
    assert len(a) == len(b)
    mp = {elt:[] for elt in np.unique(b)}
    for i in range(len(b)):
        mp[b[i]].append(a[i])
    for elt in mp:
        mp[elt] = np.unique(mp[elt]).tolist()
    return mp


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 midi2lists.py <midi_file_original> <midi_file_bp>")
        sys.exit(1)

    midi_file_1 = sys.argv[1]
    midi_file_2 = sys.argv[2]
    
    notes_with_start_time = get_notes_with_start_time(midi_file_1)
    notes_with_start_time_2 = get_notes_with_start_time(midi_file_2)

    # print("Notes with start times:", filter_tuples(notes_with_start_time))
    # print("Notes with start times:", filter_tuples(notes_with_start_time_2))
    local_maxima = find_local_maximums(filter_tuples(notes_with_start_time))
    local_maxima_2 = find_local_maximums(filter_tuples(notes_with_start_time_2))
    # print(f"Local maximums in the list: {local_maxima}")
    # print(f"Local maximums in the list: {local_maxima_2}")
    a, b = match_pitches_and_return_times(local_maxima, local_maxima_2)
    print(f"--------------------First List------------------------")
    print(f"List 1: {a}")
    print(f"--------------------Second List------------------------")
    print(f"List 2: {b}")
    print(f"--------------------Map------------------------")
    print(check_point_mapping(a, b))
    print(f"------------------------------------------------------------")

