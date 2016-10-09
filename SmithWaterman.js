var BLOSUM62 = require('./BLOSUM62');
var fs = require('fs');

var LINE_LENGTH = 60;
var NUM_P_VALUE_TRIALS = 1000;

var args = process.argv.slice(2);
var s_file = args[0];
var t_file = args[1];

var s_contents = fs.readFileSync(s_file, 'utf8');
s_contents = s_contents.split("\n")
var s_name = s_contents[0].split("|")[1];
var s = s_contents.slice(1, s_contents.length).join("");

var t_contents = fs.readFileSync(t_file, 'utf8');
t_contents = t_contents.split("\n")
var t_name = t_contents[0].split("|")[1];
var t = t_contents.slice(1, t_contents.length).join("");

var score_matrix = generateScoreMatrix(s, t);
traceBack(s, t, score_matrix);

var score = getMaxScore(score_matrix);

console.log("Alignment Score: " + score);

var p_value = getPValue(s, t, score);

console.log("p_value: " + p_value);

function generateScoreMatrix(s, t) {
    // make matrix var here
    var score_matrix = [];

    // Read sequences from files into s and t
    for (var s_index = 0; s_index < s.length; s_index++) {
        for (var t_index = 0; t_index < t.length; t_index++) {
            if (s_index == 0 || t_index == 0) {
                if (!score_matrix[s_index]) {
                    score_matrix[s_index] = [];
                }
                score_matrix[s_index][t_index] = 0;
                continue;
            }

            var diagonal_score = score_matrix[s_index - 1][t_index - 1];
            var vertical_score = score_matrix[s_index - 1][t_index];
            var horizontal_score = score_matrix[s_index][t_index - 1];

            var match_s_t = BLOSUM62.getAlignmentScore(s[s_index], t[t_index]);
            var match_s_dash = BLOSUM62.getAlignmentScore(s[s_index], "-");
            var match_dash_t = BLOSUM62.getAlignmentScore("-", t[t_index]);

            var score = Math.max(
                diagonal_score + match_s_t,
                vertical_score + match_s_dash,
                horizontal_score + match_dash_t
            );

            score_matrix[s_index][t_index] = score;
        }
    }

    return score_matrix;
}

function traceBack(s, t, score_matrix) {
    var max_score = 0;
    var max_score_s_index = 0;
    var max_score_t_index = 0;

    for (var s_index = 0; s_index < score_matrix.length; s_index++) {
        for (var t_index = 0; t_index < score_matrix[0].length; t_index++) {
            if (score_matrix[s_index][t_index] > max_score) {
                max_score = score_matrix[s_index][t_index];
                max_score_s_index = s_index;
                max_score_t_index = t_index;
            }
        }
    }

    var optimal_s_alignment = "";
    var optimal_t_alignment = "";

    var current_s_index = max_score_s_index;
    var current_t_index = max_score_t_index;
    var current_score = max_score;

    while (current_score != 0) {
        var diagonal_score = score_matrix[current_s_index - 1][current_t_index - 1];
        var vertical_score = score_matrix[current_s_index - 1][current_t_index];
        var horizontal_score = score_matrix[current_s_index][current_t_index - 1];

        var match_s_t = BLOSUM62.getAlignmentScore(s[current_s_index], t[current_t_index]);
        var match_s_dash = BLOSUM62.getAlignmentScore(s[current_s_index], "-");
        var match_dash_t = BLOSUM62.getAlignmentScore("-", t[current_t_index]);

        if (diagonal_score + match_s_t == current_score)
        {
            optimal_s_alignment = s[current_s_index] + optimal_s_alignment;
            optimal_t_alignment = t[current_t_index] + optimal_t_alignment;
            current_s_index--;
            current_t_index--;
        }
        else if (vertical_score + match_s_dash == current_score)
        {
            optimal_s_alignment = s[current_s_index] + optimal_s_alignment;
            optimal_t_alignment = "-" + optimal_t_alignment;
            current_s_index--;
        }
        else
        {
            optimal_s_alignment = "-" + optimal_s_alignment;
            optimal_t_alignment = t[current_t_index] + optimal_t_alignment;
            current_t_index--;
        }
        current_score = score_matrix[current_s_index][current_t_index];
    }

    for (var i = 0; i < Math.ceil(optimal_s_alignment.length / LINE_LENGTH); i++) {
        var end_index = Math.min((i*LINE_LENGTH) + LINE_LENGTH, optimal_s_alignment.length);

        var optimal_s_sub = optimal_s_alignment.substring(i * LINE_LENGTH, end_index);
        console.log(s_name + ":\t" + current_s_index + "\t" + optimal_s_sub);

        current_s_index += optimal_s_sub.replace("-", "").length;

        var comparison = ""
        for (var j = LINE_LENGTH * i; j < end_index; j++) {
            if (optimal_s_alignment[j] == optimal_t_alignment[j]) {
                comparison += optimal_s_alignment[j];
            } else if (BLOSUM62.getAlignmentScore(optimal_s_alignment[j], optimal_t_alignment[j]) > 0) {
                comparison += "+";
            } else {
                comparison += " ";
            }
        }
        console.log("\t\t" + comparison);


        var optimal_t_sub = optimal_t_alignment.substring(i*LINE_LENGTH, end_index);
        console.log(t_name + ":\t" + current_t_index + "\t" + optimal_t_sub);

        current_t_index += optimal_t_sub.replace("_", "").length;

        console.log();
    }
}

function getPValue(s, t, score) {
    var num_greater = 0;

    for (var i = 0; i < NUM_P_VALUE_TRIALS; i++) {
        var permuted_s = permute(s);
        var score_matrix = generateScoreMatrix(permuted_s, t);
        var max_score = getMaxScore(score_matrix);
        /*if (i % 10 == 0) {
            console.log(max_score);
        }*/

        if (max_score >= score) {
            num_greater++;
        }
    }

    var p_value = (num_greater / NUM_P_VALUE_TRIALS);

    if (p_value == 0) {
        return (1 / NUM_P_VALUE_TRIALS).toExponential();
    }

    return p_value.toExponential();
}

function getMaxScore(score_matrix) {
    var max_score = 0;

    for (var s_index = 0; s_index < score_matrix.length; s_index++) {
        for (var t_index = 0; t_index < score_matrix[0].length; t_index++) {
            if (score_matrix[s_index][t_index] > max_score) {
                max_score = score_matrix[s_index][t_index];
                max_score_s_index = s_index;
                max_score_t_index = t_index;
            }
        }
    }
    return max_score;
}

function permute(sequence) {
    var sequence_array = sequence.split("");

    for (var i = sequence_array.length - 1; i > 0; i--) {
        // Get a random index
        var j = Math.round(Math.random() * i);

        // Swap letters at indices i and j
        var temp = sequence_array[i];
        sequence_array[i] = sequence_array[j];
        sequence_array[j] = temp;
    }


    return sequence_array.join("");
}
