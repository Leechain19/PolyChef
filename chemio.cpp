//
// Created by AnthonyZhang on 2025/1/12.
//

#include "chemio.h"

std::variant<long, double, std::string> chemio::convertPyObject(PyObject* obj) {
    if (PyLong_Check(obj)) {
        return PyLong_AS_LONG(obj);
    }
    if (PyFloat_Check(obj)) {
        return PyFloat_AS_DOUBLE(obj);
    }
    if (PyUnicode_Check(obj)) {
        return std::string(PyUnicode_AsUTF8(obj));
    }
    throw std::runtime_error("Unsupported type in tuple");
}

std::vector<std::variant<long, double, std::string>> chemio::extract_tuple(PyObject* tuple) {
    std::vector<std::variant<long, double, std::string>> result;
    if (!PyTuple_Check(tuple)) {
        throw std::runtime_error("Expected a tuple");
    }

    Py_ssize_t size = PyTuple_Size(tuple);
    for (Py_ssize_t i = 0; i < size; ++i) {
        PyObject* item = PyTuple_GetItem(tuple, i);
        result.push_back(convertPyObject(item));
    }
    return result;
}

std::vector<std::vector<std::variant<long, double, std::string>>> chemio::extract_list(PyObject* list) {
    std::vector<std::vector<std::variant<long, double, std::string>>> result;
    if (!PyList_Check(list)) {
        throw std::runtime_error("Expected a list");
    }

    Py_ssize_t size = PyList_Size(list);
    for (Py_ssize_t i = 0; i < size; ++i) {
        PyObject* py_tuple = PyList_GetItem(list, i);
        if (!PyTuple_Check(py_tuple)) {
            throw std::runtime_error("List item is not a tuple");
        }
        result.push_back(extract_tuple(py_tuple));
    }
    return result;
}

std::shared_ptr<Graph> chemio::PyInfoConvertToGraph(
        const std::vector<std::vector<std::variant<long, double, std::string>>> &atoms_vec,
        const std::vector<std::vector<std::variant<long, double, std::string>>> &edges_vec,
        const std::vector<std::vector<std::variant<long, double, std::string>>> &polys_vec,
        const std::vector<std::vector<std::variant<long, double, std::string>>> &ring_edges_vec
        ) {
    auto g = std::make_shared<Graph>();

    // atoms
    for (const auto& vec : atoms_vec) {
        std::string name = std::get<std::string>(vec[0]);
        float x = static_cast<float>(std::get<double>(vec[1]));
        float y = static_cast<float>(std::get<double>(vec[2]));
        float z = static_cast<float>(std::get<double>(vec[3]));
        bool ar = static_cast<bool>(std::get<long>(vec[4]));
        g->addAtom(std::make_shared<Atom>(name, x, y, z), 0, ar);
    }

    // edges

    for (const auto& vec : edges_vec) {
        int x = static_cast<int>(std::get<long>(vec[0]));
        int y = static_cast<int>(std::get<long>(vec[1]));
        std::string type = std::get<std::string>(vec[2]);
        g->addEdge(x, y, type);
        g->addEdge(y, x, type);
    }

    // polys
    for (const auto& vec : polys_vec) {
        int who = static_cast<int>(std::get<long>(vec[0]));
        float x = static_cast<float>(std::get<double>(vec[1]));
        float y = static_cast<float>(std::get<double>(vec[2]));
        float z = static_cast<float>(std::get<double>(vec[3]));
        g->addPoly(x, y, z, who);
    }

    // ring
    for (const auto& vec : ring_edges_vec) {
        int x = static_cast<int>(std::get<long>(vec[0]));
        int y = static_cast<int>(std::get<long>(vec[1]));
        g->addRingEdge(x, y);
    }

    g->calMainChain();

    return g;
}

std::shared_ptr<Graph> chemio::buildGraphFromPSmiles(const std::string& psmiles) {
    if (psmiles.empty()) {
        throw exception::InvalidParameterException("The psmiles string is empty");
    }
    std::string smi;
    for (const auto& c : psmiles) {
        if (c == '*') {
            smi += "13C";
        }
        else smi += c;
    }

    std::cout << "psmiles = " << psmiles << std::endl;
    std::cout << "smi = " << smi << std::endl;

    // Py
    Py_Initialize();
    if (!Py_IsInitialized()) {
        std::cerr << "Python initialization failed" << std::endl;
        return nullptr;
    }

    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append('./../script')");

    PyObject* module = PyImport_ImportModule("rdkit_helper");

    if (!module) {
        PyErr_Print();
        std::cerr << "Failed to load module 'rdkit_helper.py'" << std::endl;
        Py_Finalize();
        return nullptr;
    }

    PyObject* func = PyObject_GetAttrString(module, "read_smi");
    if (!func || !PyCallable_Check(func)) {
        PyErr_Print();
        std::cerr << "Failed to load function 'read_smi'" << std::endl;
        Py_DECREF(module);
        Py_Finalize();
        return nullptr;
    }

    PyObject* args = PyTuple_New(1);
    PyTuple_SET_ITEM(args, 0, Py_BuildValue("s", smi.c_str()));

    PyObject* result = PyObject_CallObject(func, args);
    if (!result) {
        PyErr_Print();
        std::cerr << "Function call failed" << std::endl;
        Py_DECREF(func);
        Py_DECREF(module);
        Py_Finalize();
        return nullptr;
    }

    // 检查返回值是否为元组
    if (!PyTuple_Check(result)) {
        std::cerr << "Function did not return a tuple" << std::endl;
        Py_DECREF(result);
        Py_DECREF(func);
        Py_DECREF(module);
        Py_Finalize();
        return nullptr;
    }

    // 提取返回的 4 个列表
    PyObject* atoms_cache = PyTuple_GetItem(result, 0);
    PyObject* edges_cache = PyTuple_GetItem(result, 1);
    PyObject* polys_cache = PyTuple_GetItem(result, 2);
    PyObject* ring_edges_cache = PyTuple_GetItem(result, 3);

    // 转换为 C++ 数据结构
    std::vector<std::vector<std::variant<long, double, std::string>>> atoms_vec = chemio::extract_list(atoms_cache);
    std::vector<std::vector<std::variant<long, double, std::string>>> edges_vec = chemio::extract_list(edges_cache);
    std::vector<std::vector<std::variant<long, double, std::string>>> polys_vec = chemio::extract_list(polys_cache);
    std::vector<std::vector<std::variant<long, double, std::string>>> ring_edges_vec = chemio::extract_list(ring_edges_cache);

    // 清理资源
    Py_DECREF(args);
    Py_DECREF(result);
    Py_DECREF(func);
    Py_DECREF(module);

    return PyInfoConvertToGraph(atoms_vec, edges_vec, polys_vec, ring_edges_vec);
}

std::shared_ptr<Graph> chemio::buildGraphFromMol2(const std::string& path) {
    if (!std::filesystem::exists(path)) {
        throw exception::IllegalStringException("buildGraphFromMol2: Illegal Mol2 Path");
    }

    // Py
    Py_Initialize();
    if (!Py_IsInitialized()) {
        std::cerr << "Python initialization failed" << std::endl;
        return nullptr;
    }

    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append('./../script')");

    PyObject* module = PyImport_ImportModule("rdkit_helper");

    if (!module) {
        PyErr_Print();
        std::cerr << "Failed to load module 'rdkit_helper'" << std::endl;
        Py_Finalize();
        return nullptr;
    }

    PyObject* func = PyObject_GetAttrString(module, "read_mol2");
    if (!func || !PyCallable_Check(func)) {
        PyErr_Print();
        std::cerr << "Failed to load function 'read_mol2'" << std::endl;
        Py_DECREF(module);
        Py_Finalize();
        return nullptr;
    }

    PyObject* args = PyTuple_New(1);
    PyTuple_SET_ITEM(args, 0, Py_BuildValue("s", path.c_str()));

    PyObject* result = PyObject_CallObject(func, args);
    if (!result) {
        PyErr_Print();
        std::cerr << "Function call failed" << std::endl;
        Py_DECREF(func);
        Py_DECREF(module);
        Py_Finalize();
        return nullptr;
    }

    // 检查返回值是否为元组
    if (!PyTuple_Check(result)) {
        std::cerr << "Function did not return a tuple" << std::endl;
        Py_DECREF(result);
        Py_DECREF(func);
        Py_DECREF(module);
        Py_Finalize();
        return nullptr;
    }

    // 提取返回的 4 个列表
    PyObject* atoms_cache = PyTuple_GetItem(result, 0);
    PyObject* edges_cache = PyTuple_GetItem(result, 1);
    PyObject* polys_cache = PyTuple_GetItem(result, 2);
    PyObject* ring_edges_cache = PyTuple_GetItem(result, 3);

    // 转换为 C++ 数据结构
    std::vector<std::vector<std::variant<long, double, std::string>>> atoms_vec = chemio::extract_list(atoms_cache);
    std::vector<std::vector<std::variant<long, double, std::string>>> edges_vec = chemio::extract_list(edges_cache);
    std::vector<std::vector<std::variant<long, double, std::string>>> polys_vec = chemio::extract_list(polys_cache);
    std::vector<std::vector<std::variant<long, double, std::string>>> ring_edges_vec = chemio::extract_list(ring_edges_cache);

    // 清理资源
    Py_DECREF(args);
    Py_DECREF(result);
    Py_DECREF(func);
    Py_DECREF(module);

    return PyInfoConvertToGraph(atoms_vec, edges_vec, polys_vec, ring_edges_vec);
}

std::string chemio::getAtomType(const std::string& name, int bond_num, bool ar) {
    /*
    >> C.3 sp3 carbon
    >> C.2 sp2 carbon
    >> C.1 sp carbon
    >> C.ar aromatic carbon
    >> C.cat 碳正离子
     */
    if (ar) {
        return name + ".ar";
    }
    if (name == "H")
        return name;
    if (name == "C") {
        if (bond_num <= 4 && bond_num >= 2)
            return "C." + std::to_string(bond_num - 1);
        return "C.cat";
    }
    /*
    >> N.4 sp3 nitrogen 季铵盐
    >> N.3 sp3 nitrogen
    >> N.2 sp2 nitrogen
    >. N.1 sp nitrogen
    >> N.ar aromatic nitrogen 芳香氮
    >> N.am amide nitrogen 酰胺氮
     */

    if (name == "N") {
        if (bond_num > 0 && bond_num <= 4)
            return "N." + std::to_string(bond_num);
        return "N.am";
    }

    /*
    >> O.3  sp3 oxygen
    >> O.2  sp2 oxygen
     */

    if (name == "O") {
        if (bond_num == 2)
            return "O.3";
        return "O.2";
    }

    return name;
}

void chemio::writeMol2File(const std::string& file_name, const std::string& adj_file_name, const std::shared_ptr<Graph>& g, const std::string& file_info) {
    std::ofstream outFile(file_name);
    if (!outFile.is_open()) {
        std::cerr << "Error to open" << file_name << std::endl;
        return;
    }
    outFile << "###" << std::endl;
    std::time_t now = std::time(nullptr);
    char* time_str = std::ctime(&now);
    outFile << "### Created by PolyChef-cpp17 on " << time_str;
    outFile << "###" << std::endl;
    outFile << std::endl;

    outFile << "@<TRIPOS>MOLECULE" << std::endl;
    outFile << file_info << std::endl;
    outFile << g->size() << ' ' << g->getEdgesSize() << ' ' << 1 << ' ' << 0 << ' ' << 0 << std::endl;

    outFile << "SMALL" << std::endl;
    outFile << "USER_CHARGES" << std::endl;
    outFile << std::endl << std::endl;
    outFile << "@<TRIPOS>ATOM" << std::endl;

    for (int i = 0; i < g->size(); i ++) {
        const auto& ptr = g->getAtom(i);
        outFile << std::setw(7) << std::right << (i + 1) << ' ' << std::setw(8) << std::left << ptr->getSymbol() << std::setw(10) << std::fixed << std::setprecision(4) << std::right
        << ptr->getx() << std::setw(10) << std::fixed << std::setprecision(4) << std::right << ptr->gety() << std::setw(10) << std::fixed << std::setprecision(4) << std::right
        << ptr->getz() << ' ' << std::setw(5) << std::left << chemio::getAtomType(ptr->getSymbol(), (int)g->getEdge(i).size(), g->isAr(i)) << std::setw(6) << std::right
        << g->getMonomer(i) << std::setw(9) << std::right << 0 << std::setw(10) << std::fixed <<  std::setprecision(4) << std::right << 0.0f << '\n';
    }

    outFile << "@<TRIPOS>BOND" << std::endl;
    int edge_idx = 1;

    for (int i = 0; i < (int)g->getEdgesVec().size(); i ++) {
        const auto& vec = g->getEdge(i);
        for (const auto& ptr : vec) {
            int to = ptr->getTo();
            if (i < to) {
                const std::string& tp = ptr->getType();
                outFile << std::setw(6) << std::right << (edge_idx ++) << std::setw(6) << std::right << i + 1 << std::setw(6) << std::right
                << to + 1 << std::setw(7) << std::right << tp << '\n';
            }
        }
    }
    outFile.close();

    std::ofstream outFileAdj(adj_file_name);
    for (int i = 0; i < g->size(); i ++) {
        std::string s = std::to_string(i + 1);
        for (auto& edgePtr : g->getEdge(i)) {
            int to = edgePtr->getTo();
            s.push_back(' ');
            s += g->getAtom(to)->getSymbol();
        }
        outFileAdj << s << std::endl;
    }
    outFileAdj.close();

    std::cout << "The file named " << file_name << " has been finished" << std::endl;
    std::cout << "The file named " << adj_file_name << " has been finished" << std::endl;
}

void chemio::writeLoss2File(const std::string& loss_file_name, const std::unique_ptr<std::vector<std::pair<double, double>>>& loss_vec_ptr) {
    std::ofstream outFile(loss_file_name);
    if (!outFile.is_open()) {
        std::cerr << "Error to open" << loss_file_name << std::endl;
        return;
    }

    for (const auto& p : *loss_vec_ptr) {
        outFile << p.first << ' ' << p.second << '\n';
    }
    outFile << std::endl;
    outFile.close();

    std::cout << "The file named " << loss_file_name << " has been finished" << std::endl;
}
