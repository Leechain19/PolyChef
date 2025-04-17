//
// Created by AnthonyZhang on 2025/1/12.
//

#include "chemio.h"

Index2CodePrinter::Index2CodePrinter(std::string header) : header_(std::move(header)) {}

std::string Index2CodePrinter::index2code(int idx) {
    assert(idx > 0);
    std::string ans;
    while (idx > 0) {
        auto u = idx % 26;
        if (u == 0) u = 26;
        ans.push_back(static_cast<char>('A' + u - 1));
        idx = (idx - u) / 26;
    }
    std::reverse(ans.begin(), ans.end());
    return ans;
}

std::string Index2CodePrinter::get(int memo_type) {
    while (memo_type >= (int)this->string_memo.size()) {
        this->string_memo.emplace_back();
    }
    if (this->string_memo[memo_type].empty()) {
        this->string_memo[memo_type] = header_ + index2code(memo_type);
    }
    return this->string_memo[memo_type];
}

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
        g->addAtom(std::make_shared<Atom>(name, x, y, z), 0, 0, ar);
    }

    // edges

    for (const auto& vec : edges_vec) {
        int x = static_cast<int>(std::get<long>(vec[0]));
        int y = static_cast<int>(std::get<long>(vec[1]));
        std::string bond_string = std::get<std::string>(vec[2]);
        Bond_type bond_type = getBondTypeByString(bond_string);
        g->addEdge(x, y, bond_type);
        g->addEdge(y, x, bond_type);
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

std::shared_ptr<CrossLinker> chemio::PyInfoConvertToCrossLinker(
        const std::vector<std::vector<std::variant<long, double, std::string>>> &atoms_vec,
        const std::vector<std::vector<std::variant<long, double, std::string>>> &edges_vec,
        const std::vector<std::vector<std::variant<long, double, std::string>>> &polys_vec,
        const std::vector<std::vector<std::variant<long, double, std::string>>> &ring_edges_vec
        ) {
    auto g = std::make_shared<CrossLinker>();

    // atoms
    for (const auto& vec : atoms_vec) {
        std::string name = std::get<std::string>(vec[0]);
        float x = static_cast<float>(std::get<double>(vec[1]));
        float y = static_cast<float>(std::get<double>(vec[2]));
        float z = static_cast<float>(std::get<double>(vec[3]));
        bool ar = static_cast<bool>(std::get<long>(vec[4]));
        g->addAtom(std::make_shared<Atom>(name, x, y, z), ar);
    }

    // edges

    for (const auto& vec : edges_vec) {
        int x = static_cast<int>(std::get<long>(vec[0]));
        int y = static_cast<int>(std::get<long>(vec[1]));
        std::string bond_string = std::get<std::string>(vec[2]);
        Bond_type bond_type = chemio::getBondTypeByString(bond_string);

        g->addEdge(x, y, bond_type);
        g->addEdge(y, x, bond_type);
    }

    // polys
    for (const auto& vec : polys_vec) {
        int who = static_cast<int>(std::get<long>(vec[0]));
        float x = static_cast<float>(std::get<double>(vec[1]));
        float y = static_cast<float>(std::get<double>(vec[2]));
        float z = static_cast<float>(std::get<double>(vec[3]));
        g->addPoly(x, y, z, who);
    }

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
    const char* pythonScriptPath = PYTHON_SCRIPT_PATH;
    std::string simple_string = std::string("sys.path.append('") + std::string(pythonScriptPath) + std::string("')");
    PyRun_SimpleString(simple_string.c_str());

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
    const char* pythonScriptPath = PYTHON_SCRIPT_PATH;
    std::string simple_string = std::string("sys.path.append('") + std::string(pythonScriptPath) + std::string("')");
    PyRun_SimpleString(simple_string.c_str());

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

std::shared_ptr<CrossLinker> chemio::buildCrossLinkerFromPSmiles(const std::string& psmiles) {
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
    const char* pythonScriptPath = PYTHON_SCRIPT_PATH;
    std::string simple_string = std::string("sys.path.append('") + std::string(pythonScriptPath) + std::string("')");
    PyRun_SimpleString(simple_string.c_str());

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

    return PyInfoConvertToCrossLinker(atoms_vec, edges_vec, polys_vec, ring_edges_vec);
}

std::shared_ptr<CrossLinker> chemio::buildCrossLinkerFromMol2(const std::string& path) {
    if (!std::filesystem::exists(path)) {
        throw exception::IllegalStringException("buildCrossLinkerFromMol2: Illegal Mol2 Path");
    }
    // Py
    Py_Initialize();
    if (!Py_IsInitialized()) {
        std::cerr << "Python initialization failed" << std::endl;
        return nullptr;
    }

    PyRun_SimpleString("import sys");
    const char* pythonScriptPath = PYTHON_SCRIPT_PATH;
    std::string simple_string = std::string("sys.path.append('") + std::string(pythonScriptPath) + std::string("')");
    PyRun_SimpleString(simple_string.c_str());

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

    return PyInfoConvertToCrossLinker(atoms_vec, edges_vec, polys_vec, ring_edges_vec);
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

std::string chemio::getBondSymbolByBondType(Bond_type bond_type) {
    if (bond_type == Bond_type::SINGLE_BOND) {
        return "1";
    }
    if (bond_type == Bond_type::DOUBLE_BOND) {
        return "2";
    }
    if (bond_type == Bond_type::TRIPLE_BOND) {
        return "3";
    }
    return "ar";
}

Bond_type chemio::getBondTypeByString(const std::string& s) {
    if (s == "1") return Bond_type::SINGLE_BOND;
    if (s == "2") return Bond_type::DOUBLE_BOND;
    if (s == "3") return Bond_type::TRIPLE_BOND;
    return Bond_type::AROMATIC_BOND;
}

void chemio::writeMol2File(const std::string& file_name, const std::string& adj_file_name, const std::shared_ptr<Graph>& g, const std::string& file_info) {
    std::ofstream outFile(file_name);
    if (!outFile.is_open()) {
        std::cerr << "Error to open" << file_name << std::endl;
        return;
    }
    Index2CodePrinter mono_printer(std::string("M_"));

    outFile << "###" << std::endl;
    std::time_t now = std::time(nullptr);
    char* time_str = std::ctime(&now);
    outFile << "### Created by PolyChef-cpp17 on " << time_str;
    outFile << "###" << std::endl;
    outFile << std::endl;

    outFile << "@<TRIPOS>MOLECULE" << std::endl;
    outFile << file_info << std::endl;
    outFile << g->size() << ' ' << g->getEdgesSum() << ' ' << 1 << ' ' << 0 << ' ' << 0 << std::endl;

    outFile << "SMALL" << std::endl;
    outFile << "USER_CHARGES" << std::endl;
    outFile << std::endl << std::endl;
    outFile << "@<TRIPOS>ATOM" << std::endl;

    for (int i = 0; i < g->size(); i ++) {
        const auto& ptr = g->getAtom(i);
        outFile << std::setw(7) << std::right << (i + 1) << ' ' << std::setw(8) << std::left << ptr->getSymbol() << std::setw(10) << std::fixed << std::setprecision(4) << std::right
        << ptr->getx() << std::setw(10) << std::fixed << std::setprecision(4) << std::right << ptr->gety() << std::setw(10) << std::fixed << std::setprecision(4) << std::right
        << ptr->getz() << ' ' << std::setw(5) << std::left << chemio::getAtomType(ptr->getSymbol(), (int)g->getEdge(i).size(), g->isAr(i)) << std::setw(6) << std::right
        << g->getMonomer(i) << std::setw(9) << std::right << mono_printer.get(g->getMonomerType(i)) << std::setw(10) << std::fixed <<  std::setprecision(4) << std::right << 0.0f << '\n';
    }

    outFile << "@<TRIPOS>BOND" << std::endl;
    int edge_idx = 1;

    for (int i = 0; i < (int)g->getEdgesVec().size(); i ++) {
        const auto& vec = g->getEdge(i);
        for (const auto& ptr : vec) {
            int to = ptr->getTo();
            if (i < to) {
                auto tp = ptr->getType();
                outFile << std::setw(6) << std::right << (edge_idx ++) << std::setw(6) << std::right << i + 1 << std::setw(6) << std::right
                << to + 1 << std::setw(7) << std::right << getBondSymbolByBondType(tp) << '\n';
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

void chemio::writeMol2File(const std::string& file_name, const std::string& adj_file_name, const std::unique_ptr<CrosslinkingSystem>& cls, const std::string& file_info) {
    std::ofstream outFile(file_name);
    if (!outFile.is_open()) {
        std::cerr << "Error to open " << file_name << std::endl;
        return;
    }

    Index2CodePrinter mono_printer(std::string("M_"));
    Index2CodePrinter cross_printer(std::string("C_"));

    outFile << "###\n";
    std::time_t now = std::time(nullptr);
    char* time_str = std::ctime(&now);
    outFile << "### Created by PolyChef-cpp17 on " << time_str;
    outFile << "###\n";
    outFile << "\n";

    outFile << "@<TRIPOS>MOLECULE\n";
    outFile << file_info << std::endl;
    int atom_sum = cls->getAtomSize();
    int edge_sum = cls->getEdgeSize();

    outFile << atom_sum << ' ' << edge_sum << ' ' << 1 << ' ' << 0 << ' ' << 0 << std::endl;

    outFile << "SMALL\n";
    outFile << "USER_CHARGES\n";
    outFile << "\n\n";
    outFile << "@<TRIPOS>ATOM\n" << std::flush;


    // crosslinker first
    std::vector<std::vector<int>> chain_index_table(cls->getCrosslinkerNetworkNumber());
    std::vector<std::vector<int>> crosslink_index_table(cls->getCrosslinkerNumber());
    std::vector<std::string> atom_env(atom_sum);

    int mono_index = 0, atom_index = 0, edge_index = 0;

    int cl_sum_atom = 0;

    for (int gid = 0; gid < cls->getCrosslinkerNumber(); gid ++) {
        const auto& cl = cls->getCrosslinkGraph(gid);
        cl_sum_atom += cl->size();
        for (int i = 0; i < cl->size(); i ++) {
            const auto atom_ptr = cl->getAtom(i);

            outFile << std::setw(7) << std::right << (atom_index + 1) << ' ' << std::setw(8) << std::left << atom_ptr->getSymbol() << std::setw(10) << std::fixed << std::setprecision(4) << std::right
            << atom_ptr->getx() << std::setw(10) << std::fixed << std::setprecision(4) << std::right << atom_ptr->gety() << std::setw(10) << std::fixed << std::setprecision(4) << std::right
            << atom_ptr->getz() << ' ' << std::setw(5) << std::left << chemio::getAtomType(atom_ptr->getSymbol(), (int)cl->getEdge(i).size(), cl->isAr(i)) << std::setw(6) << std::right
            << mono_index << std::setw(9) << std::right << cross_printer.get(cl->getMonomerType()) << std::setw(10) << std::fixed <<  std::setprecision(4) << std::right << 0.0f << '\n';

            crosslink_index_table[gid].emplace_back(atom_index);
            atom_index ++;
        }
        mono_index ++;
    }

    // chain graphs
    for (int gid = 0; gid < cls->getCrosslinkerNetworkNumber(); gid ++) {
        auto chain = cls->getChainGraph(gid);
        int next_mono = mono_index;
        for (int i = 0; i < chain->size(); i ++) {
            const auto atom_ptr = chain->getAtom(i);

            outFile << std::setw(7) << std::right << (atom_index + 1) << ' ' << std::setw(8) << std::left << atom_ptr->getSymbol() << std::setw(10) << std::fixed << std::setprecision(4) << std::right
            << atom_ptr->getx() << std::setw(10) << std::fixed << std::setprecision(4) << std::right << atom_ptr->gety() << std::setw(10) << std::fixed << std::setprecision(4) << std::right
            << atom_ptr->getz() << ' ' << std::setw(5) << std::left << chemio::getAtomType(atom_ptr->getSymbol(), (int)chain->getEdge(i).size(), chain->isAr(i)) << std::setw(6) << std::right
            << (chain->getMonomer(i) + mono_index) << std::setw(9) << std::right << mono_printer.get(chain->getMonomerType(i)) << std::setw(10) << std::fixed <<  std::setprecision(4) << std::right << 0.0f << '\n';

            next_mono = std::max(next_mono, chain->getMonomer(i) + mono_index);

            chain_index_table[gid].emplace_back(atom_index);
            atom_index ++;
        }
        mono_index = next_mono + 1;
    }

    // bond
    outFile << "@<TRIPOS>BOND" << std::endl;
    for (int gid = 0; gid < cls->getCrosslinkerNumber(); gid ++) {
        const auto& cl = cls->getCrosslinkGraph(gid);
        for (int i = 0; i < cl->size(); i ++) {
            const auto& vec = cl->getEdge(i);
            auto atom_id1 = crosslink_index_table[gid][i];

            for (const auto& e : vec) {
                auto to = e->getTo();
                auto atom_id2 = crosslink_index_table[gid][to];
                atom_env[atom_id1] += cl->getAtom(to)->getSymbol() + " ";

                if (atom_id1 < atom_id2) {
                    auto tp = e->getType();
                    outFile << std::setw(6) << std::right << (++ edge_index) << std::setw(6) << std::right << atom_id1 + 1 << std::setw(6) << std::right
                        << atom_id2 + 1 << std::setw(7) << std::right << getBondSymbolByBondType(tp) << '\n';
                }
            }
        }
    }

    for (int gid = 0; gid < cls->getCrosslinkerNetworkNumber(); gid ++) {
        const auto& g = cls->getChainGraph(gid);
        for (int i = 0; i < g->size(); i ++) {
            const auto& vec = g->getEdge(i);
            auto atom_id1 = chain_index_table[gid][i];

            for (const auto& e : vec) {
                auto to = e->getTo();
                auto atom_id2 = chain_index_table[gid][to];
                atom_env[atom_id1] += g->getAtom(to)->getSymbol() + " ";

                if (atom_id1 < atom_id2) {
                    auto tp = e->getType();
                    outFile << std::setw(6) << std::right << (++ edge_index) << std::setw(6) << std::right << atom_id1 + 1 << std::setw(6) << std::right <<
                        atom_id2 + 1 << std::setw(7) << std::right << chemio::getBondSymbolByBondType(tp) << '\n';
                }
            }
        }
    }

    // linker
    for (int i = 0; i < cls->getCrosslinkerNetworkNumber(); i ++) {
        const auto& e = cls->getCrosslinkerNetworkInIdx(i);
        {
            int cross_id = e[0], cross_poly_id = e[1];
            auto atom_id1 = crosslink_index_table[cross_id][cls->getCrosslinkGraph(cross_id)->getPoly(cross_poly_id)->getNeigh()];
            auto chain_id = i;
            auto atom_id2 = chain_index_table[chain_id][cls->getChainGraph(i)->polyFront()->getNeigh()];

            atom_env[atom_id1] += cls->getChainGraph(i)->getAtom(cls->getChainGraph(i)->polyFront()->getNeigh())->getSymbol() + " ";
            atom_env[atom_id2] += cls->getCrosslinkGraph(cross_id)->getAtom(cls->getCrosslinkGraph(cross_id)->getPoly(cross_poly_id)->getNeigh())->getSymbol() + " ";

            outFile << std::setw(6) << std::right << (++ edge_index) << std::setw(6) << std::right << atom_id1 + 1 << std::setw(6) << std::right
            << atom_id2 + 1 << std::setw(7) << std::right << "1" << '\n';
        }
        {
            int cross_id = e[2], cross_poly_id = e[3];
            auto atom_id1 = crosslink_index_table[cross_id][cls->getCrosslinkGraph(cross_id)->getPoly(cross_poly_id)->getNeigh()];
            auto chain_id = i;
            auto atom_id2 = chain_index_table[chain_id][cls->getChainGraph(i)->polyBack()->getNeigh()];

            atom_env[atom_id1] += cls->getChainGraph(i)->getAtom(cls->getChainGraph(i)->polyBack()->getNeigh())->getSymbol() + " ";
            atom_env[atom_id2] += cls->getCrosslinkGraph(cross_id)->getAtom(cls->getCrosslinkGraph(cross_id)->getPoly(cross_poly_id)->getNeigh())->getSymbol() + " ";

            outFile << std::setw(6) << std::right << (++ edge_index) << std::setw(6) << std::right << atom_id1 + 1 << std::setw(6) << std::right
            << atom_id2 + 1 << std::setw(7) << std::right << "1" << '\n';
        }
    }

    std::cout << "atom_idx: " << atom_index << " atom_sum: " << atom_sum << " edge_idx: " << edge_index << " edge_sum: " << edge_sum << std::endl;
    outFile.close();

    assert(atom_index == atom_sum && edge_index == edge_sum);

    std::ofstream adjFile(adj_file_name);
    if (!adjFile.is_open()) {
        std::cerr << "Error to open " << adj_file_name << std::endl;
        outFile.close();
        return;
    }

    for (int idx = 0; idx < (int)atom_env.size(); idx ++) {
        adjFile << idx + 1 << ' ' << atom_env[idx] << '\n';
    }
    adjFile.close();

    std::cout << "The file named " << file_name << " has been finished" << std::endl;
    std::cout << "The file named " << adj_file_name << " has been finished" << std::endl;
}