var BLOSUM62 = require('./BLOSUM62');
var fs = require('fs');

var LINE_LENGTH = 60;

//var s = "ARNDCQEKMFPSTWYVBZX"; // first sequence
//var t = "ARDCQEGHILMFAAAASTWYX"; // second sequence
var args = process.argv.slice(2);
var s_file = args[0];
var t_file = args[1];

var s = fs.readFileSync(s_file, 'utf8');
s = s.split("\n").slice(1, s.length).join("");
var t = fs.readFileSync(t_file, 'utf8');
t = t.split("\n").slice(1, t.length).join("");

var result = generateScoreMatrix(s, t);
traceBack(s, t, result);

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
        console.log("s:\t" + current_s_index + "\t" + optimal_s_sub);

        current_s_index += optimal_s_sub.replace("-", "").length;

        var comparison = ""
        for (var j = LINE_LENGTH * i; j < end_index; j++) {
            if (optimal_s_alignment[j] != "-" && optimal_s_alignment[j] == optimal_t_alignment[j]) {
                comparison += optimal_s_alignment[j];
            } else if (BLOSUM62.getAlignmentScore(optimal_s_alignment[j], optimal_t_alignment[j]) > 0) {
                comparison += "+";
            } else {
                comparison += " ";
            }
        }
        console.log("\t\t" + comparison);


        var optimal_t_sub = optimal_t_alignment.substring(i*LINE_LENGTH, end_index);
        console.log("t:\t" + current_t_index + "\t" + optimal_t_sub);

        current_t_index += optimal_t_sub.replace("_", "").length;

        console.log();
    }
}
