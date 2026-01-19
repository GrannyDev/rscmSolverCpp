//
// Snapshot persistence for recomputing costs without re-solving.
//

#include "SnapshotIO.h"
#include "../solver/Solver.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <boost/dynamic_bitset.hpp>

namespace {

std::string BitsetToString(const boost::dynamic_bitset<>& bs) {
    std::string out;
    boost::to_string(bs, out);
    return out;
}

boost::dynamic_bitset<> BitsetFromString(const std::string& s) {
    return boost::dynamic_bitset<>(s);
}

template <typename T>
std::string Join(const std::vector<T>& v, char sep = ',') {
    std::ostringstream oss;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) oss << sep;
        oss << v[i];
    }
    return oss.str();
}

template <typename T>
std::vector<T> Split(const std::string& s) {
    std::vector<T> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (tok.empty()) continue;
        std::stringstream conv(tok);
        T val;
        conv >> val;
        out.push_back(val);
    }
    return out;
}

} // namespace

void WriteSnapshot(const Solver& solver, const RSCM& solutionNode, const std::string& path, bool overwrite)
{
    if (std::filesystem::exists(path) && !overwrite) {
        throw std::runtime_error("Snapshot already exists: " + path);
    }
    std::ofstream ofs(path, std::ios::trunc);
    if (!ofs) throw std::runtime_error("Unable to open snapshot file: " + path);

    ofs << "layout:" << Join(solver.layout) << "\n";
    ofs << "targets:" << Join(solver.targets) << "\n";
    ofs << "nb_input_bits:" << solver.nbInputBits << "\n";
    ofs << "lut_width:" << solver.lutWidth_ << "\n";
    ofs << "max_coef:" << solver.maxCoef << "\n";
    ofs << "min_coef:" << solver.minCoef << "\n";
    ofs << "nb_bits_per_scm:" << solver.nbBitsPerSCM << "\n";
    ofs << "nb_adders:" << solver.nbAdders << "\n";
    ofs << "nb_variables:" << solver.nbPossibleVariables << "\n";
    ofs << "nb_muxes:" << solver.varDefs.size() << "\n";
    ofs << "is_symmetric:" << (solver.isSymmetric_ ? 1 : 0) << "\n";

    ofs << "scm_indexes:" << Join(solutionNode.scmIndexes) << "\n";
    ofs << "rscm_set:" << BitsetToString(solutionNode.rscm.set) << "\n";
    ofs << "rscm_max:" << Join(solutionNode.rscm.maxOutputValue) << "\n";
    ofs << "rscm_min:" << Join(solutionNode.rscm.minOutputValue) << "\n";
    ofs << "rscm_coeff_tz:" << Join(solutionNode.rscm.coefficientTrailingZeros) << "\n";
    ofs << "rscm_is_left_minus:" << Join(solutionNode.rscm.isLeftMinus, ',') << "\n";
    ofs << "rscm_is_right_minus:" << Join(solutionNode.rscm.isRightMinus, ',') << "\n";
    ofs << "min_shift_savings:" << Join(solutionNode.minShiftSavings) << "\n";
    ofs << "variable_bit_widths:" << Join(solutionNode.variableBitWidths) << "\n";
    ofs << "is_plus_minus:" << Join(solutionNode.isPlusMinus, ',') << "\n";

    // Store selected SCM DAGs (one per target)
    ofs << "selected_scms:" << solver.targets.size() << "\n";
    for (size_t i = 0; i < solver.targets.size(); ++i) {
        const unsigned idx = solutionNode.scmIndexes[i];
        const auto& dag = solver.scmDesigns[i].second[idx];
        ofs << "scm_target:" << solver.scmDesigns[i].first << "\n";
        ofs << "scm_set:" << BitsetToString(dag.set) << "\n";
        ofs << "scm_max:" << Join(dag.maxOutputValue) << "\n";
        ofs << "scm_min:" << Join(dag.minOutputValue) << "\n";
        ofs << "scm_coeff_tz:" << Join(dag.coefficientTrailingZeros) << "\n";
        ofs << "scm_is_left_minus:" << Join(dag.isLeftMinus, ',') << "\n";
        ofs << "scm_is_right_minus:" << Join(dag.isRightMinus, ',') << "\n";
    }
}

std::optional<SnapshotData> ReadSnapshot(const std::string& path)
{
    std::ifstream ifs(path);
    if (!ifs) return std::nullopt;
    SnapshotData snap;
    std::string line;
    size_t expectedScms = 0;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        auto colon = line.find(':');
        if (colon == std::string::npos) continue;
        const std::string key = line.substr(0, colon);
        const std::string val = line.substr(colon + 1);
        if (key == "layout") snap.layout = Split<int>(val);
        else if (key == "targets") snap.targets = Split<int>(val);
        else if (key == "nb_input_bits") snap.nbInputBits = static_cast<size_t>(std::stoul(val));
        else if (key == "lut_width") snap.lutWidth = static_cast<size_t>(std::stoul(val));
        else if (key == "max_coef") snap.maxCoef = std::stoi(val);
        else if (key == "min_coef") snap.minCoef = std::stoi(val);
        else if (key == "nb_bits_per_scm") snap.nbBitsPerSCM = static_cast<size_t>(std::stoul(val));
        else if (key == "nb_adders") snap.nbAdders = static_cast<size_t>(std::stoul(val));
        else if (key == "nb_variables") snap.nbVariables = static_cast<size_t>(std::stoul(val));
        else if (key == "nb_muxes") snap.nbMuxes = static_cast<size_t>(std::stoul(val));
        else if (key == "is_symmetric") snap.isSymmetric = (std::stoi(val) != 0);
        else if (key == "scm_indexes") snap.scmIndexes = Split<unsigned int>(val);
        else if (key == "rscm_set") snap.rscmSet = BitsetFromString(val);
        else if (key == "rscm_max") snap.rscmMaxOutput = Split<int>(val);
        else if (key == "rscm_min") snap.rscmMinOutput = Split<int>(val);
        else if (key == "rscm_coeff_tz") snap.rscmCoeffTZ = Split<unsigned int>(val);
        else if (key == "rscm_is_left_minus") {
            auto tmp = Split<int>(val);
            snap.rscmIsLeftMinus.assign(tmp.begin(), tmp.end());
        }
        else if (key == "rscm_is_right_minus") {
            auto tmp = Split<int>(val);
            snap.rscmIsRightMinus.assign(tmp.begin(), tmp.end());
        }
        else if (key == "rscm_is_minus") {
            auto tmp = Split<int>(val);
            snap.rscmIsRightMinus.assign(tmp.begin(), tmp.end());
        }
        else if (key == "min_shift_savings") snap.minShiftSavings = Split<unsigned int>(val);
        else if (key == "variable_bit_widths") snap.variableBitWidths = Split<unsigned int>(val);
        else if (key == "is_plus_minus") {
            auto tmp = Split<int>(val);
            snap.isPlusMinus.assign(tmp.begin(), tmp.end());
        }
        else if (key == "selected_scms") expectedScms = static_cast<size_t>(std::stoul(val));
        else if (key == "scm_target") {
            DAG dag(static_cast<unsigned int>(snap.nbBitsPerSCM), static_cast<unsigned int>(snap.nbAdders), static_cast<unsigned int>(snap.nbVariables));
            int target = std::stoi(val);
            std::string setLine, maxLine, minLine, coeffLine, minusLine;
            if (!std::getline(ifs, setLine) || !std::getline(ifs, maxLine) || !std::getline(ifs, minLine) || !std::getline(ifs, coeffLine) || !std::getline(ifs, minusLine)) break;
            dag.set = BitsetFromString(setLine.substr(setLine.find(':') + 1));
            dag.maxOutputValue = Split<int>(maxLine.substr(maxLine.find(':') + 1));
            dag.minOutputValue = Split<int>(minLine.substr(minLine.find(':') + 1));
            dag.coefficientTrailingZeros = Split<unsigned int>(coeffLine.substr(coeffLine.find(':') + 1));
            const auto minusKey = minusLine.substr(0, minusLine.find(':'));
            if (minusKey == "scm_is_left_minus") {
                auto tmpLeft = Split<int>(minusLine.substr(minusLine.find(':') + 1));
                dag.isLeftMinus.assign(tmpLeft.begin(), tmpLeft.end());
                std::string rightLine;
                if (!std::getline(ifs, rightLine)) break;
                auto tmpRight = Split<int>(rightLine.substr(rightLine.find(':') + 1));
                dag.isRightMinus.assign(tmpRight.begin(), tmpRight.end());
            } else {
                auto tmpRight = Split<int>(minusLine.substr(minusLine.find(':') + 1));
                dag.isRightMinus.assign(tmpRight.begin(), tmpRight.end());
                dag.isLeftMinus.assign(tmpRight.size(), false);
            }
            snap.selectedScms.emplace_back(std::move(dag));
            (void)target; // stored in same order as targets
        }
    }
    if (expectedScms && snap.selectedScms.size() != expectedScms) {
        return std::nullopt;
    }
    if (snap.rscmIsLeftMinus.empty() && !snap.rscmIsRightMinus.empty()) {
        snap.rscmIsLeftMinus.assign(snap.rscmIsRightMinus.size(), false);
    }
    if (snap.rscmIsRightMinus.empty() && !snap.rscmIsLeftMinus.empty()) {
        snap.rscmIsRightMinus.assign(snap.rscmIsLeftMinus.size(), false);
    }
    return snap;
}
