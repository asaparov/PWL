#ifndef SET_GRAPH_H_
#define SET_GRAPH_H_

#include <core/array.h>
#include <core/map.h>
#include <stdint.h>
#include <set>

using namespace core;


/* forward declarations */

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer> struct set_reasoning;

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
bool find_largest_disjoint_clique_with_set(
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>&,
		unsigned int, unsigned int*&, unsigned int&, unsigned int&, int);
template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
bool find_largest_disjoint_clique_with_set(
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>&,
		unsigned int, unsigned int, unsigned int*&, unsigned int&, unsigned int&, int);

template<typename T, typename Stream>
bool print(const hash_set<T>& set, Stream&& out) {
	if (!print('{', out)) return false;
	bool first = true;
	for (const T& element : set) {
		if (first) {
			first = false;
		} else if (!print(", ", out)) {
			return false;
		}
		if (!print(element, out)) return false;
	}
	return print('}', out);
}

struct intensional_set_vertex {
	array<unsigned int> parents;
	array<unsigned int> children;

	static inline bool clone(const intensional_set_vertex& src, intensional_set_vertex& dst) {
		if (!array_init(dst.parents, src.parents.capacity)) {
			return false;
		} else if (!array_init(dst.children, src.children.capacity)) {
			core::free(dst.parents);
			return false;
		}
		for (unsigned int parent : src.parents)
			dst.parents[dst.parents.length++] = parent;
		for (unsigned int child : src.children)
			dst.children[dst.children.length++] = child;
		return true;
	}

	static inline void free(intensional_set_vertex& vertex) {
		core::free(vertex.parents);
		core::free(vertex.children);
	}
};

inline bool init(intensional_set_vertex& vertex) {
	if (!array_init(vertex.parents, 4)) {
		return false;
	} else if (!array_init(vertex.children, 4)) {
		core::free(vertex.parents);
		return false;
	}
	return true;
}

template<typename Stream>
bool read(intensional_set_vertex& vertex, Stream& in)
{
	if (!read(vertex.parents, in)) {
		return false;
	} else if (!read(vertex.children, in)) {
		free(vertex.parents);
		return false;
	}
	return true;
}

template<typename Stream>
bool write(const intensional_set_vertex& vertex, Stream& out)
{
	return write(vertex.parents, out)
		&& write(vertex.children, out);
}

template<typename ProofCalculus>
struct extensional_set_vertex
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;

	array_map<unsigned int, array<Proof*>> parents;
	array_map<unsigned int, array<Proof*>> children;

	/* NOTE: we can't properly implement a `clone` function because this struct
	   only owns the memory of the proofs in `children` and not those in
	   `parents`, so to correctly clone an extensional graph, we first need to
	   clone all of the proofs in `children` across all vertices, and then
	   clone the proofs in `parents`. */
	static inline bool clone_except_parents(
			const extensional_set_vertex<ProofCalculus>& src,
			extensional_set_vertex<ProofCalculus>& dst,
			array_map<const Proof*, Proof*>& proof_map,
			hash_map<const Formula*, Formula*>& formula_map)
	{
		if (!array_map_init(dst.parents, src.parents.capacity)) {
			return false;
		} else if (!array_map_init(dst.children, src.children.capacity)) {
			core::free(dst.parents);
			return false;
		}
		for (const auto& entry : src.children) {
			dst.children.keys[dst.children.size] = entry.key;
			array<Proof*>& dst_proofs = dst.children.values[dst.children.size];
			if (!array_init(dst_proofs, entry.value.capacity)) {
				core::free(dst);
				return false;
			}
			dst.children.size++;

			for (Proof* proof : entry.value) {
				if (!Proof::clone(proof, dst_proofs[dst_proofs.length], proof_map, formula_map)) {
					core::free(dst);
					return false;
				}
				dst_proofs.length++;
			}
		}
		return true;
	}

	static inline bool clone_only_parents(
			const extensional_set_vertex<ProofCalculus>& src,
			extensional_set_vertex<ProofCalculus>& dst,
			array_map<const Proof*, Proof*>& proof_map)
	{
		for (const auto& entry : src.parents) {
			dst.parents.keys[dst.parents.size] = entry.key;
			array<Proof*>& dst_proofs = dst.parents.values[dst.parents.size];
			if (!array_init(dst_proofs, entry.value.capacity))
				return false;
			dst.parents.size++;

			for (Proof* proof : entry.value) {
				unsigned int index = proof_map.index_of(proof);
#if !defined(NDEBUG)
				if (index == proof_map.size)
					fprintf(stderr, "extensional_set_vertex.clone_only_parents WARNING: Given proof does not exist in `proof_map`.\n");
#endif
				dst_proofs[dst_proofs.length++] = proof_map.values[index];
			}
		}
		return true;
	}

	static inline void free(extensional_set_vertex<ProofCalculus>& vertex) {
		for (auto entry : vertex.parents)
			core::free(entry.value);
		for (auto entry : vertex.children) {
			for (Proof* proof : entry.value) {
				core::free(*proof); if (proof->reference_count == 0) core::free(proof);
			}
			core::free(entry.value);
		}
		core::free(vertex.parents);
		core::free(vertex.children);
	}
};

template<typename ProofCalculus>
inline bool init(extensional_set_vertex<ProofCalculus>& vertex) {
	if (!array_map_init(vertex.parents, 4)) {
		return false;
	} else if (!array_map_init(vertex.children, 4)) {
		core::free(vertex.parents);
		return false;
	}
	return true;
}

template<typename ProofCalculus>
inline bool get_proof_map(
		const extensional_set_vertex<ProofCalculus>& vertex,
		hash_map<const typename ProofCalculus::Proof*, unsigned int>& proof_map,
		hash_map<const typename ProofCalculus::Language*, unsigned int>& formula_map)
{
	typedef typename ProofCalculus::Proof Proof;

	for (const auto& entry : vertex.parents) {
		for (const Proof* proof : entry.value) {
			if (!get_proof_map(proof, proof_map, formula_map))
				return false;
		}
	} for (const auto& entry : vertex.children) {
		for (const Proof* proof : entry.value) {
			if (!get_proof_map(proof, proof_map, formula_map))
				return false;
		}
	}
	return true;
}

template<typename ProofCalculus, typename Stream>
bool read(
		extensional_set_vertex<ProofCalculus>& vertex,
		Stream& in,
		typename ProofCalculus::Proof** proofs)
{
	decltype(vertex.parents.size) parent_count;
	decltype(vertex.children.size) child_count;

	if (!read(parent_count, in)
	 || !read(child_count, in))
		return false;

	if (!array_map_init(vertex.parents, ((size_t) 1) << (core::log2(parent_count == 0 ? 1 : parent_count) + 1))) {
		return false;
	} else if (!array_map_init(vertex.children, ((size_t) 1) << (core::log2(child_count == 0 ? 1 : child_count) + 1))) {
		free(vertex.parents);
		return false;
	}

	for (unsigned int i = 0; i < parent_count; i++) {
		size_t length;
		if (!read(vertex.parents.keys[i], in)
		 || !read(length, in)
		 || !array_init(vertex.parents.values[i], ((size_t) 1) << (core::log2(length == 0 ? 1 : length) + 1)))
		{
			free(vertex);
			return false;
		}
		vertex.parents.size++;
		for (size_t j = 0; j < length; j++) {
			unsigned int index;
			if (!read(index, in)) {
				free(vertex);
				return false;
			}
			vertex.parents.values[i][j] = proofs[index];
			vertex.parents.values[i].length++;
		}
	} for (unsigned int i = 0; i < child_count; i++) {
		size_t length;
		if (!read(vertex.children.keys[i], in)
		 || !read(length, in)
		 || !array_init(vertex.children.values[i], ((size_t) 1) << (core::log2(length == 0 ? 1 : length) + 1)))
		{
			free(vertex);
			return false;
		}
		vertex.children.size++;
		for (size_t j = 0; j < length; j++) {
			unsigned int index;
			if (!read(index, in)) {
				free(vertex);
				return false;
			}
			vertex.children.values[i][j] = proofs[index];
			proofs[index]->reference_count++;
			vertex.children.values[i].length++;
		}
	}
	return true;
}

template<typename ProofCalculus, typename Stream>
bool write(
		const extensional_set_vertex<ProofCalculus>& vertex, Stream& out,
		const hash_map<const typename ProofCalculus::Proof*, unsigned int>& proof_map)
{
	typedef typename ProofCalculus::Proof Proof;

	if (!write(vertex.parents.size, out)
	 || !write(vertex.children.size, out))
		return false;
	for (const auto& entry : vertex.parents) {
		if (!write(entry.key, out)
		 || !write(entry.value.length, out))
			return false;
		for (const Proof* proof : entry.value)
			if (!write(proof_map.get(proof), out)) return false;
	} for (const auto& entry : vertex.children) {
		if (!write(entry.key, out)
		 || !write(entry.value.length, out))
			return false;
		for (const Proof* proof : entry.value)
			if (!write(proof_map.get(proof), out)) return false;
	}
	return true;
}

struct intensional_set_graph {
	intensional_set_vertex* vertices;

	intensional_set_graph(unsigned int initial_capacity) {
		vertices = (intensional_set_vertex*) malloc(sizeof(intensional_set_vertex) * initial_capacity);
		if (vertices == NULL) exit(EXIT_FAILURE);
	}

	~intensional_set_graph() {
		core::free(vertices);
	}

	static inline void free(intensional_set_graph& graph) {
		core::free(graph.vertices);
	}

	inline bool resize(unsigned int new_capacity) {
		return core::resize(vertices, new_capacity);
	}

	inline bool new_set(unsigned int vertex_id) {
		return init(vertices[vertex_id]);
	}

	template<bool CleanupEdges>
	inline void free_set(unsigned int vertex_id) {
		if (CleanupEdges) {
			for (unsigned int parent : vertices[vertex_id].parents)
				remove_edge(parent, vertex_id);
			for (unsigned int child : vertices[vertex_id].children)
				remove_edge(vertex_id, child);
		}
		core::free(vertices[vertex_id]);
	}

	inline bool add_edge(unsigned int parent, unsigned int child) {
		if (!vertices[parent].children.add(child)) {
			return false;
		} else if (!vertices[child].parents.add(parent)) {
			vertices[parent].children.length--;
			return false;
		}
		return true;
	}

	inline void remove_edge(unsigned int parent, unsigned int child) {
		unsigned int index = vertices[parent].children.index_of(child);
#if !defined(NDEBUG)
		if (index == vertices[parent].children.length)
			fprintf(stderr, "intensional_set_graph.remove_edge WARNING: The set %u is not in `vertices[%u].children`.\n", child, parent);
#endif
		vertices[parent].children.remove(index);
		index = vertices[child].parents.index_of(parent);
#if !defined(NDEBUG)
		if (index == vertices[child].parents.length)
			fprintf(stderr, "intensional_set_graph.remove_edge WARNING: The set %u is not in `vertices[%u].children`.\n", parent, child);
#endif
		vertices[child].parents.remove(index);
	}
};

template<typename ProofCalculus>
struct extensional_set_graph
{
	typedef typename ProofCalculus::Proof Proof;

	extensional_set_vertex<ProofCalculus>* vertices;

	extensional_set_graph(unsigned int initial_capacity) {
		vertices = (extensional_set_vertex<ProofCalculus>*) malloc(sizeof(extensional_set_vertex<ProofCalculus>) * initial_capacity);
		if (vertices == NULL) exit(EXIT_FAILURE);
	}

	~extensional_set_graph() {
		core::free(vertices);
	}

	static inline void free(extensional_set_graph<ProofCalculus>& graph) {
		core::free(graph.vertices);
	}

	inline bool resize(unsigned int new_capacity) {
		return core::resize(vertices, new_capacity);
	}

	inline bool new_set(unsigned int vertex_id) {
		return init(vertices[vertex_id]);
	}

	template<bool CleanupEdges>
	inline void free_set(unsigned int vertex_id) {
		if (CleanupEdges) {
			for (const auto& entry : vertices[vertex_id].parents) {
				core::free(entry.value);
			} for (const auto& entry : vertices[vertex_id].children) {
				for (Proof* proof : entry.value)
					core::free(*proof);
				core::free(entry.value);
			}
		}
		core::free(vertices[vertex_id]);
	}

	inline bool add_edge(unsigned int parent, unsigned int child, Proof* axiom)
	{
		if (!vertices[parent].children.ensure_capacity(vertices[parent].children.size + 1)
		 || !vertices[child].parents.ensure_capacity(vertices[child].parents.size + 1))
			return false;

		unsigned int index = vertices[parent].children.index_of(child);
		array<Proof*>& parent_axioms = vertices[parent].children.values[index];
		if (index == vertices[parent].children.size) {
			if (!array_init(parent_axioms, 4)) return false;
			vertices[parent].children.keys[index] = child;
			vertices[parent].children.size++;
		} else if (!parent_axioms.ensure_capacity(parent_axioms.length + 1)) {
			return false;
		}

		index = vertices[child].parents.index_of(parent);
		array<Proof*>& child_axioms = vertices[child].parents.values[index];
		if (index == vertices[child].parents.size) {
			if (!array_init(child_axioms, 4)) return false;
			vertices[child].parents.keys[index] = parent;
			vertices[child].parents.size++;
		} else if (!child_axioms.ensure_capacity(child_axioms.length + 1)) {
			return false;
		}

#if !defined(NDEBUG)
		for (Proof* existing_axiom : parent_axioms) {
			if (existing_axiom == axiom || *existing_axiom == *axiom) {
				fprintf(stderr, "extensional_edge.add_edge WARNING: The given edge already exists.\n");
				return false;
			}
		} for (Proof* existing_axiom : child_axioms) {
			if (existing_axiom == axiom || *existing_axiom == *axiom) {
				fprintf(stderr, "extensional_edge.add_edge WARNING: The given edge already exists.\n");
				return false;
			}
		}
#endif

		parent_axioms[parent_axioms.length++] = axiom;
		child_axioms[child_axioms.length++] = axiom;
		axiom->reference_count++;
		return true;
	}

	template<typename Formula>
	inline Proof* get_existing_edge(unsigned int parent, unsigned int child,
			Formula* parent_formula, Formula* child_formula) const
	{
		typedef typename Formula::Type FormulaType;

		array<Proof*>& axioms = vertices[parent].children.get(child);
		for (Proof* axiom : axioms) {
			Formula* formula = axiom->formula->quantifier.operand;
			while (formula->type == FormulaType::FOR_ALL)
				formula = formula->quantifier.operand;
			if ((formula->binary.left == child_formula || *formula->binary.left == *child_formula)
			 && (formula->binary.right == parent_formula || *formula->binary.right == *parent_formula))
			{
				return axiom;
			}
		}
		return nullptr;
	}

	template<typename Formula>
	inline Proof* get_edge(unsigned int parent, unsigned int child,
			Formula* parent_formula, Formula* child_formula,
			unsigned int arity, bool& new_edge)
	{
		typedef typename Formula::Type FormulaType;

		if (!vertices[parent].children.ensure_capacity(vertices[parent].children.size + 1)
		 || !vertices[child].parents.ensure_capacity(vertices[child].parents.size + 1))
			return nullptr;

		unsigned int index = vertices[parent].children.index_of(child);
		array<Proof*>& parent_axioms = vertices[parent].children.values[index];
		if (index == vertices[parent].children.size) {
			if (!array_init(parent_axioms, 4)) return nullptr;
			vertices[parent].children.keys[index] = child;
			vertices[parent].children.size++;
		} else if (!parent_axioms.ensure_capacity(parent_axioms.length + 1)) {
			return nullptr;
		}

		index = vertices[child].parents.index_of(parent);
		array<Proof*>& child_axioms = vertices[child].parents.values[index];
		if (index == vertices[child].parents.size) {
			if (!array_init(child_axioms, 4)) return nullptr;
			vertices[child].parents.keys[index] = parent;
			vertices[child].parents.size++;
		} else if (!child_axioms.ensure_capacity(child_axioms.length + 1)) {
			return nullptr;
		}

		for (Proof* existing_axiom : parent_axioms) {
			Formula* formula = existing_axiom->formula->quantifier.operand;
			while (formula->type == FormulaType::FOR_ALL)
				formula = formula->quantifier.operand;
			if ((formula->binary.left == child_formula || *formula->binary.left == *child_formula)
			 && (formula->binary.right == parent_formula || *formula->binary.right == *parent_formula))
			{
				new_edge = false;
				return existing_axiom;
			}
		}

		new_edge = true;
		Formula* formula = Formula::new_for_all(arity, Formula::new_if_then(child_formula, parent_formula));
		if (formula == nullptr) return nullptr;
		for (unsigned int i = arity - 1; i > 0; i--) {
			Formula* new_formula = Formula::new_for_all(i, formula);
			if (new_formula == nullptr) {
				core::free(*formula); core::free(formula);
				return nullptr;
			}
			formula = new_formula;
		}
		Proof* new_axiom = ProofCalculus::new_axiom(formula);
		core::free(*formula); if (formula->reference_count == 0) core::free(formula);
		if (new_axiom == nullptr) return nullptr;
		child_formula->reference_count++;
		parent_formula->reference_count++;

		parent_axioms[parent_axioms.length++] = new_axiom;
		child_axioms[child_axioms.length++] = new_axiom;
		new_axiom->reference_count++;
		return new_axiom;
	}

	template<typename Formula>
	inline void remove_edge(unsigned int parent, unsigned int child,
			Formula* parent_formula, Formula* child_formula)
	{
		typedef typename Formula::Type FormulaType;

		unsigned int index = vertices[child].parents.index_of(parent);
#if !defined(NDEBUG)
		if (index == vertices[child].parents.size)
			fprintf(stderr, "extensional_set_graph.remove_edge WARNING: The set %u is not in `vertices[%u].children`.\n", child, parent);
#endif
		array<Proof*>& child_axioms = vertices[child].parents.values[index];
		for (unsigned int i = 0; i < child_axioms.length; i++) {
			Proof* axiom = child_axioms[i];
			Formula* formula = axiom->formula->quantifier.operand;
			while (formula->type == FormulaType::FOR_ALL)
				formula = formula->quantifier.operand;
			if ((formula->binary.left == child_formula || *formula->binary.left == *child_formula)
			 && (formula->binary.right == parent_formula || *formula->binary.right == *parent_formula))
			{
				child_axioms.remove(i);
				break;
			}
		}
		if (child_axioms.length == 0) {
			core::free(child_axioms);
			vertices[child].parents.remove_at(index);
		}

		index = vertices[parent].children.index_of(child);
#if !defined(NDEBUG)
		if (index == vertices[parent].children.size)
			fprintf(stderr, "extensional_set_graph.remove_edge WARNING: The set %u is not in `vertices[%u].children`.\n", parent, child);
#endif
		array<Proof*>& parent_axioms = vertices[parent].children.values[index];
		for (unsigned int i = 0; i < parent_axioms.length; i++) {
			Proof* axiom = parent_axioms[i];
			Formula* formula = axiom->formula->quantifier.operand;
			while (formula->type == FormulaType::FOR_ALL)
				formula = formula->quantifier.operand;
			if ((formula->binary.left == child_formula || *formula->binary.left == *child_formula)
			 && (formula->binary.right == parent_formula || *formula->binary.right == *parent_formula))
			{
				core::free(*axiom);
				if (axiom->reference_count == 0)
					core::free(axiom);
				parent_axioms.remove(i);
				break;
			}
		}
		if (parent_axioms.length == 0) {
			core::free(parent_axioms);
			vertices[parent].children.remove_at(index);
		}
	}
};

bool get_ancestors(
		const intensional_set_graph& graph, unsigned int vertex,
		array_map<unsigned int, unsigned int>& ancestors)
{
	if (!ancestors.ensure_capacity(ancestors.size + 1)) return false;
	size_t index = ancestors.index_of(vertex);
	if (index < ancestors.size) {
		ancestors.values[index]++;
		return true;
	}
	ancestors.keys[index] = vertex;
	ancestors.values[index] = 0;
	ancestors.size++;

	array<unsigned int> stack(8);
	stack.add(vertex);
	while (stack.length > 0) {
		unsigned int v = stack.pop();
		for (unsigned int parent : graph.vertices[v].parents) {
			if (!ancestors.ensure_capacity(ancestors.size + 1)) return false;
			size_t index = ancestors.index_of(parent);
			if (index < ancestors.size) {
				ancestors.values[index]++;
			} else {
				ancestors.keys[index] = parent;
				ancestors.values[index] = 1;
				ancestors.size++;
				if (!stack.add(parent)) return false;
			}
		}
	}
	return true;
}

bool get_descendants(
		const intensional_set_graph& graph, unsigned int vertex,
		array_map<unsigned int, unsigned int>& descendants)
{
	if (!descendants.ensure_capacity(descendants.size + 1)) return false;
	size_t index = descendants.index_of(vertex);
	if (index < descendants.size) {
		descendants.values[index]++;
		return true;
	}
	descendants.keys[index] = vertex;
	descendants.values[index] = 0;
	descendants.size++;

	array<unsigned int> stack(8);
	stack.add(vertex);
	while (stack.length > 0) {
		unsigned int v = stack.pop();
		for (unsigned int child : graph.vertices[v].children) {
			if (!descendants.ensure_capacity(descendants.size + 1)) return false;
			size_t index = descendants.index_of(child);
			if (index < descendants.size) {
				descendants.values[index]++;
			} else {
				descendants.keys[index] = child;
				descendants.values[index] = 1;
				descendants.size++;
				if (!stack.add(child)) return false;
			}
		}
	}
	return true;
}

typedef uint_fast8_t tuple_element_type_specifier;
enum class tuple_element_type : tuple_element_type_specifier {
	CONSTANT,
	NUMBER,
	STRING
};

struct tuple_element {
	tuple_element_type type;
	union {
		unsigned int constant;
		hol_number number;
		string str;
	};

	static inline unsigned int hash(const tuple_element& key) {
		/* TODO: precompute these statically */
		unsigned int type_hash = default_hash(key.type);
		switch (key.type) {
		case tuple_element_type::CONSTANT:
			return type_hash ^ default_hash(key.constant);
		case tuple_element_type::NUMBER:
			return type_hash ^ default_hash(key.number);
		case tuple_element_type::STRING:
			return type_hash ^ string::hash(key.str);
		}
		fprintf(stderr, "tuple_element.hash ERROR: Unrecognized `tuple_element_type`.\n");
		exit(EXIT_FAILURE);
	}

	static inline void move(const tuple_element& src, tuple_element& dst) {
		dst.type = src.type;
		switch (src.type) {
		case tuple_element_type::CONSTANT:
			dst.constant = src.constant; return;
		case tuple_element_type::NUMBER:
			dst.number = src.number; return;
		case tuple_element_type::STRING:
			core::move(src.str, dst.str); return;
		}
		fprintf(stderr, "tuple_element.move ERROR: Unrecognized `tuple_element_type`.\n");
		exit(EXIT_FAILURE);
	}

	static inline void free(tuple_element& element) {
		if (element.type == tuple_element_type::STRING)
			core::free(element.str);
	}

private:
	inline bool init_helper(const tuple_element& src) {
		type = src.type;
		switch (type) {
		case tuple_element_type::CONSTANT:
			constant = src.constant; return true;
		case tuple_element_type::NUMBER:
			number = src.number; return true;
		case tuple_element_type::STRING:
			return init(str, src.str);
		}
		fprintf(stderr, "tuple_element.init_helper ERROR: Unrecognized `tuple_element_type`.\n");
		return false;
	}

	friend bool init(tuple_element&, const tuple_element&);
};

inline bool init(tuple_element& element, const tuple_element& src) {
	return element.init_helper(src);
}

inline bool operator == (const tuple_element& first, const tuple_element& second) {
	if (first.type != second.type)
		return false;
	switch (first.type) {
	case tuple_element_type::CONSTANT:
		return first.constant == second.constant;
	case tuple_element_type::NUMBER:
		return first.number == second.number;
	case tuple_element_type::STRING:
		return first.str == second.str;
	}
	fprintf(stderr, "operator == ERROR: Unrecognized `tuple_element_type`.\n");
	return false;
}

inline bool operator != (const tuple_element& first, const tuple_element& second) {
	if (first.type != second.type)
		return true;
	switch (first.type) {
	case tuple_element_type::CONSTANT:
		return first.constant != second.constant;
	case tuple_element_type::NUMBER:
		return first.number != second.number;
	case tuple_element_type::STRING:
		return first.str != second.str;
	}
	fprintf(stderr, "operator != ERROR: Unrecognized `tuple_element_type`.\n");
	return false;
}

inline bool operator < (const tuple_element& first, const tuple_element& second) {
	if (first.type < second.type) return true;
	else if (first.type > second.type) return false;
	switch (first.type) {
	case tuple_element_type::CONSTANT:
		return first.constant < second.constant;
	case tuple_element_type::NUMBER:
		return first.number < second.number;
	case tuple_element_type::STRING:
		return first.str < second.str;
	}
	fprintf(stderr, "operator < ERROR: Unrecognized `tuple_element_type`.\n");
	return false;
}


template<typename Stream, typename... Printer>
inline bool print(const tuple_element& element, Stream& out, Printer&&... constant_printer) {
	switch (element.type) {
	case tuple_element_type::CONSTANT:
		return print(element.constant, out, std::forward<Printer>(constant_printer)...);
	case tuple_element_type::NUMBER:
		return print(element.number, out);
	case tuple_element_type::STRING:
		return print('"', out) && print(element.str, out) && print('"', out);
	}
	fprintf(stderr, "print ERROR: Unrecognized `tuple_element_type`.\n");
	return false;
}

template<typename Stream>
bool read(tuple_element& tup, Stream& in)
{
	tuple_element_type_specifier type;
	if (!read(type, in))
		return false;
	tup.type = (tuple_element_type) type;
	switch (tup.type) {
	case tuple_element_type::CONSTANT:
		return read(tup.constant, in);
	case tuple_element_type::NUMBER:
		return read(tup.number, in);
	case tuple_element_type::STRING:
		return read(tup.str, in);
	}
	fprintf(stderr, "read ERROR: Unrecognized `tuple_element_type`.\n");
	return false;
}

template<typename Stream>
bool write(const tuple_element& tup, Stream& out)
{
	if (!write((tuple_element_type_specifier) tup.type, out))
		return false;
	switch (tup.type) {
	case tuple_element_type::CONSTANT:
		return write(tup.constant, out);
	case tuple_element_type::NUMBER:
		return write(tup.number, out);
	case tuple_element_type::STRING:
		return write(tup.str, out);
	}
	fprintf(stderr, "write ERROR: Unrecognized `tuple_element_type`.\n");
	return false;
}

struct tuple {
	tuple_element* elements;
	unsigned int length;

	inline bool operator = (const tuple& src) {
		return init_helper(src.elements, src.length);
	}

	inline tuple_element& operator [] (unsigned int i) {
		return elements[i];
	}

	inline const tuple_element& operator [] (unsigned int i) const {
		return elements[i];
	}

	static inline bool is_empty(const tuple& src) {
		return src.elements == nullptr;
	}

	static inline void set_empty(tuple& key) {
		key.elements = nullptr;
	}

	static inline void move(const tuple& src, tuple& dst) {
		dst.elements = src.elements;
		dst.length = src.length;
	}

	static inline void swap(tuple& first, tuple& second) {
		core::swap(first.elements, second.elements);
		core::swap(first.length, second.length);
	}

	static inline unsigned int hash(const tuple& key) {
		unsigned int value = 0;
		for (unsigned int i = 0; i < key.length; i++)
			value ^= tuple_element::hash(key.elements[i]);
		return value;
	}

	static inline void free(tuple& tup) {
		for (unsigned int i = 0; i < tup.length; i++)
			core::free(tup.elements[i]);
		core::free(tup.elements);
	}

private:
	inline bool init_helper(unsigned int src_length) {
		length = src_length;
		elements = (tuple_element*) malloc(sizeof(tuple_element) * length);
		if (elements == nullptr) {
			fprintf(stderr, "tuple.init_helper ERROR: Out of memory.\n");
			return false;
		}
		return true;
	}

	inline bool init_helper(const tuple_element* src_elements, unsigned int src_length) {
		if (!init_helper(src_length))
			return false;
		for (unsigned int i = 0; i < src_length; i++) {
			if (!init(elements[i], src_elements[i])) {
				for (unsigned int j = 0; j < i; j++)
					core::free(elements[j]);
				core::free(elements);
				return false;
			}
		}
		return true;
	}

	friend bool init(tuple&, unsigned int);
	friend bool init(tuple&, const tuple&);
	friend bool init(tuple&, const tuple_element*, unsigned int);
};

inline bool init(tuple& new_tuple, unsigned int length) {
	return new_tuple.init_helper(length);
}

inline bool init(tuple& new_tuple, const tuple& src) {
	return new_tuple.init_helper(src.elements, src.length);
}

inline bool init(tuple& new_tuple, const tuple_element* src_elements, unsigned int src_length) {
	return new_tuple.init_helper(src_elements, src_length);
}

inline bool operator == (const tuple& first, const tuple& second) {
	if (first.length != second.length)
		return false;
	/* `first` may be uninitialized */
	if (first.elements == nullptr)
		return false;
	for (unsigned int i = 0; i < first.length; i++)
		if (first.elements[i] != second.elements[i]) return false;
	return true;
}

inline bool operator != (const tuple& first, const tuple& second) {
	if (first.length != second.length)
		return true;
	/* `first` may be uninitialized */
	if (first.elements == nullptr)
		return true;
	for (unsigned int i = 0; i < first.length; i++)
		if (first.elements[i] != second.elements[i]) return true;
	return false;
}

inline bool operator < (const tuple& first, const tuple& second) {
	if (first.length < second.length) return true;
	else if (first.length > second.length) return false;
	for (unsigned int i = 0; i < first.length; i++) {
		if (first.elements[i] < second.elements[i]) return true;
		else if (second.elements[i] < first.elements[i]) return false;
	}
	return false;
}

inline bool operator >= (const tuple& first, const tuple& second) {
	return !(first < second);
}

template<typename Stream, typename... Printer>
bool print(const tuple& tup, Stream& out, Printer&&... printer) {
	if (tup.length == 0)
		return print("()", out);
	if (tup.length == 1)
		return print(tup[0], out, std::forward<Printer>(printer)...);

	if (!print('(', out) || !print(tup[0], out, std::forward<Printer>(printer)...))
		return false;
	for (unsigned int i = 1; i < tup.length; i++) {
		if (!print(", ", out) || !print(tup[i], out, std::forward<Printer>(printer)...))
			return false;
	}
	return print(')', out);
}

template<typename Stream>
bool read(tuple& tup, Stream& in)
{
	if (!read(tup.length, in))
		return false;
	tup.elements = (tuple_element*) malloc(sizeof(tuple_element) * tup.length);
	if (tup.elements == nullptr) {
		fprintf(stderr, "read ERROR: Insufficient memory for `tuple.elements`.\n");
		return false;
	}
	for (unsigned int i = 0; i < tup.length; i++) {
		if (!read(tup.elements[i], in)) {
			for (unsigned int j = 0; j < i; j++)
				free(tup.elements[j]);
			free(tup.elements);
			return false;
		}
	}
	return true;
}

template<typename Stream>
bool write(const tuple& tup, Stream& out) {
	if (!write(tup.length, out))
		return false;
	for (unsigned int i = 0; i < tup.length; i++)
		if (!write(tup.elements[i], out)) return false;
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus>
struct set_info
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename Formula::Type FormulaType;

	unsigned int arity;
	unsigned int set_size;
	array<Proof*> size_axioms;
	hash_set<unsigned int> descendants;
	array<tuple_element> elements; /* NOTE: this does not contain the elements of the descendants of this set */
	array<tuple> provable_elements;

	/* optimization: this is a list of pairs of sets that would become disjoint, if this set is empty */
	array<pair<unsigned int, unsigned int>> newly_disjoint_cache;

	inline unsigned int element_count() const {
		return elements.length / arity;
	}

	inline bool add_element(const tuple& element) {
		if (!elements.ensure_capacity(elements.length + arity))
			return false;
		for (unsigned int j = 0; j < arity; j++) {
			if (!init(elements[elements.length + j], element[j])) {
				for (unsigned int k = 0; k < j; k++)
					core::free(elements[elements.length + k]);
				return false;
			}
		}
		elements.length += arity;
		return true;
	}

	inline unsigned int index_of_element(const tuple& element) const {
		unsigned int num_elements = element_count();
		if (element.length != arity) return num_elements;
		for (unsigned int i = 0; i < num_elements; i++) {
			bool equivalent = true;
			for (unsigned int j = 0; equivalent && j < arity; j++)
				if (elements[i * arity + j] != element[j]) equivalent = false;
			if (equivalent) return i;
		}
		return num_elements;
	}

	inline void remove_element_at(unsigned int i) {
		unsigned int last = element_count() - 1;
		for (unsigned int j = 0; j < arity; j++) {
			core::free(elements[i * arity + j]);
			move(elements[last * arity + j], elements[i * arity + j]);
		}
		elements.length -= arity;
	}

	static inline bool clone(
			const set_info<BuiltInConstants, ProofCalculus>& src,
			set_info<BuiltInConstants, ProofCalculus>& dst,
			array_map<const Proof*, Proof*>& proof_map,
			hash_map<const Formula*, Formula*>& formula_map)
	{
		dst.arity = src.arity;
		dst.set_size = src.set_size;
		if (!array_init(dst.size_axioms, src.size_axioms.capacity)) {
			return false;
		} else if (!hash_set_init(dst.descendants, src.descendants.capacity)) {
			core::free(dst.size_axioms);
			return false;
		} else if (!array_init(dst.elements, src.elements.capacity)) {
			core::free(dst.size_axioms);
			core::free(dst.descendants);
			return false;
		} else if (!array_init(dst.provable_elements, src.provable_elements.capacity)) {
			core::free(dst.size_axioms);
			core::free(dst.descendants);
			core::free(dst.elements);
			return false;
		} else if (!array_init(dst.newly_disjoint_cache, src.newly_disjoint_cache.capacity)) {
			core::free(dst.provable_elements);
			core::free(dst.size_axioms);
			core::free(dst.descendants);
			core::free(dst.elements);
			return false;
		}
		for (Proof* axiom : src.size_axioms) {
			if (!Proof::clone(axiom, dst.size_axioms[dst.size_axioms.length], proof_map, formula_map)) {
				core::free(dst);
				return false;
			}
			dst.size_axioms.length++;
		} for (unsigned int descendant : src.descendants) {
			dst.descendants.add(descendant);
		} for (const tuple_element& element : src.elements) {
			if (!init(dst.elements[dst.elements.length], element)) {
				core::free(dst);
				return false;
			}
			dst.elements.length++;
		} for (const tuple& element : src.provable_elements) {
			if (!init(dst.provable_elements[dst.provable_elements.length], element)) {
				core::free(dst);
				return false;
			}
			dst.provable_elements.length++;
		} for (const pair<unsigned int, unsigned int>& p : src.newly_disjoint_cache)
			dst.newly_disjoint_cache[dst.newly_disjoint_cache.length++] = p;
		return true;
	}

	static inline void free(set_info<BuiltInConstants, ProofCalculus>& info) {
		core::free(info.newly_disjoint_cache);
		core::free(info.descendants);
		for (unsigned int j = 0; j < info.elements.length; j++)
			core::free(info.elements[j]);
		core::free(info.elements);
		for (tuple& tup : info.provable_elements)
			core::free(tup);
		core::free(info.provable_elements);
		for (Proof* axiom : info.size_axioms) {
			core::free(*axiom); if (axiom->reference_count == 0) core::free(axiom);
		}
		core::free(info.size_axioms);
		info.size_axioms.data = nullptr;
	}

	inline Formula* set_formula() {
		Formula* operand = size_axioms[0]->formula->binary.left->binary.right->quantifier.operand;
		while (operand->type == FormulaType::LAMBDA)
			operand = operand->quantifier.operand;
		return operand;
	}

	/* NOTE: this function does not check the consistency of the new size */
	inline void change_size(unsigned int new_size) {
#if !defined(NDEBUG)
		for (Proof* size_axiom : size_axioms) {
			if (size_axiom->children.length != 0 && new_size != set_size) {
				fprintf(stderr, "set_info.change_size WARNING: Attempted to change the size of a fixed set.\n");
				break;
			}
		}
#endif
		set_size = new_size;
		for (Proof* size_axiom : size_axioms)
			size_axiom->formula->binary.right->number.integer = new_size;
	}

	template<typename... Args>
	inline Proof* get_size_axiom(Formula* set_formula, Args&&... visitor) {
		if (!size_axioms.ensure_capacity(size_axioms.length + 1))
			return nullptr;
		for (Proof* size_axiom : size_axioms) {
			Formula* size_axiom_formula = size_axiom->formula->binary.left->binary.right;
			for (unsigned int i = 0; i < arity; i++)
				size_axiom_formula = size_axiom_formula->quantifier.operand;
			if (size_axiom_formula == set_formula || *size_axiom_formula == *set_formula)
				return size_axiom;
		}

		Formula* set_definition = Formula::new_lambda(arity, set_formula);
		if (set_definition == nullptr) return nullptr;
		set_formula->reference_count++;
		for (unsigned int i = arity - 1; i > 0; i--) {
			Formula* temp = Formula::new_lambda(i, set_definition);
			if (temp == NULL) {
				core::free(*set_definition); core::free(set_definition);
				return nullptr;
			}
			set_definition = temp;
		}
		Formula* size_axiom_formula = Formula::new_equals(Formula::new_atom(
				(unsigned int) BuiltInConstants::SIZE, set_definition), Formula::new_number(set_size, 0));
		if (size_axiom_formula == nullptr) {
			core::free(*set_definition); core::free(set_definition);
			return nullptr;
		}
		Proof* new_size_axiom = ProofCalculus::new_axiom(size_axiom_formula);
		core::free(*size_axiom_formula); if (size_axiom_formula->reference_count == 0) core::free(size_axiom_formula);
		if (new_size_axiom == nullptr) return nullptr;
		new_size_axiom->reference_count++;
		size_axioms[size_axioms.length++] = new_size_axiom;
		on_new_size_axiom(new_size_axiom, std::forward<Args>(visitor)...);
		return new_size_axiom;
	}

	/* NOTE: this function does not check the consistency of the new size */
	template<typename... Args>
	inline bool set_size_axiom(Proof* axiom, Args&&... visitor) {
		if (!size_axioms.ensure_capacity(size_axioms.length + 1))
			return false;
		set_size = axiom->formula->binary.right->number.integer;
		for (unsigned int i = 0; i < size_axioms.length; i++) {
			if (axiom == size_axioms[i]) {
				return true;
			} else if (*axiom->formula->binary.left->binary.right == *size_axioms[i]->formula->binary.left->binary.right) {
				core::free(*size_axioms[i]); if (size_axioms[i]->reference_count == 0) core::free(size_axioms[i]);
				size_axioms[i] = axiom;
				axiom->reference_count++;
				for (unsigned int j = 0; j < size_axioms.length; j++)
					if (i != j) size_axioms[j]->formula->binary.right->number.integer = set_size;
				return true;
			}
		}
		size_axioms[size_axioms.length++] = axiom;
		axiom->reference_count++;
		for (unsigned int j = 0; j + 1 < size_axioms.length; j++)
			size_axioms[j]->formula->binary.right->number.integer = set_size;
		on_new_size_axiom(axiom, std::forward<Args>(visitor)...);
		return true;
	}

	inline bool is_size_axiom_used_in_proof() const {
		for (Proof* size_axiom : size_axioms)
			if (size_axiom->children.length != 0) return true;
		return false;
	}
};

template<typename BuiltInConstants, typename ProofCalculus, typename... Args>
inline bool init(
		set_info<BuiltInConstants, ProofCalculus>& info,
		unsigned int arity, unsigned int set_size,
		typename ProofCalculus::Language* set_formula,
		Args&&... visitor)
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;

	info.arity = arity;
	info.set_size = set_size;
	Formula* set_definition = Formula::new_lambda(arity, set_formula);
	if (set_definition == NULL) return false;
	set_formula->reference_count++;
	for (unsigned int i = arity - 1; i > 0; i--) {
		Formula* temp = Formula::new_lambda(i, set_definition);
		if (temp == NULL) {
			free(*set_definition); free(set_definition);
			return false;
		}
		set_definition = temp;
	}
	Formula* size_axiom_formula = Formula::new_equals(Formula::new_atom(
			(unsigned int) BuiltInConstants::SIZE, set_definition), Formula::new_number(set_size, 0));
	if (size_axiom_formula == NULL) {
		free(*set_definition); free(set_definition);
		return false;
	}
	Proof* initial_size_axiom = ProofCalculus::new_axiom(size_axiom_formula);
	free(*size_axiom_formula); if (size_axiom_formula->reference_count == 0) free(size_axiom_formula);
	if (initial_size_axiom == nullptr) return false;
	initial_size_axiom->reference_count++;

	if (!array_init(info.size_axioms, 2)) {
		free(*initial_size_axiom); free(initial_size_axiom);
		return false;
	}
	info.size_axioms[info.size_axioms.length++] = initial_size_axiom;

	if (!hash_set_init(info.descendants, 16)) {
		free(*initial_size_axiom); free(initial_size_axiom); free(info.size_axioms);
		return false;
	} else if (!array_init(info.elements, 4)) {
		free(*initial_size_axiom); free(initial_size_axiom); free(info.size_axioms);
		free(info.descendants); return false;
	} else if (!array_init(info.provable_elements, 8)) {
		free(*initial_size_axiom); free(initial_size_axiom); free(info.size_axioms);
		free(info.descendants); free(info.elements); return false;
	} else if (!array_init(info.newly_disjoint_cache, 8)) {
		free(*initial_size_axiom); free(initial_size_axiom); free(info.size_axioms);
		free(info.descendants); free(info.elements); free(info.provable_elements); return false;
	}
	on_new_size_axiom(initial_size_axiom, std::forward<Args>(visitor)...);
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus>
inline bool get_proof_map(
		const set_info<BuiltInConstants, ProofCalculus>& info,
		hash_map<const typename ProofCalculus::Proof*, unsigned int>& proof_map,
		hash_map<const typename ProofCalculus::Language*, unsigned int>& formula_map)
{
	typedef typename ProofCalculus::Proof Proof;

	for (const Proof* proof : info.size_axioms)
		if (!get_proof_map(proof, proof_map, formula_map)) return false;
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Stream>
bool read(
		set_info<BuiltInConstants, ProofCalculus>& info,
		Stream& in,
		typename ProofCalculus::Proof** proofs)
{
	decltype(info.size_axioms.length) size_axiom_count;
	decltype(info.descendants.size) descendant_count;
	decltype(info.elements.length) element_count;
	decltype(info.provable_elements.length) provable_element_count;
	decltype(info.newly_disjoint_cache.length) newly_disjoint_cache_size;

	if (!read(info.arity, in)
	 || !read(info.set_size, in)
	 || !read(size_axiom_count, in)
	 || !read(descendant_count, in)
	 || !read(element_count, in)
	 || !read(provable_element_count, in)
	 || !read(newly_disjoint_cache_size, in))
		return false;

	if (!array_init(info.size_axioms, ((size_t) 1) << (core::log2(size_axiom_count == 0 ? 1 : size_axiom_count) + 1))) {
		return false;
	} else if (!hash_set_init(info.descendants, 1 << (core::log2(RESIZE_THRESHOLD_INVERSE * (descendant_count == 0 ? 1 : descendant_count)) + 1))) {
		free(info.size_axioms);
		return false;
	} else if (!array_init(info.elements, ((size_t) 1) << (core::log2(element_count == 0 ? 1 : element_count) + 1))) {
		free(info.size_axioms); free(info.descendants);
		return false;
	} else if (!array_init(info.provable_elements, ((size_t) 1) << (core::log2(provable_element_count == 0 ? 1 : provable_element_count) + 1))) {
		free(info.size_axioms); free(info.descendants);
		free(info.elements); return false;
	} else if (!array_init(info.newly_disjoint_cache, ((size_t) 1) << (core::log2(newly_disjoint_cache_size == 0 ? 1 : newly_disjoint_cache_size) + 1))) {
		free(info.size_axioms); free(info.descendants);
		free(info.elements); free(info.provable_elements);
		return false;
	}

	unsigned int index;
	for (unsigned int i = 0; i < size_axiom_count; i++) {
		if (!read(index, in)) {
			free(info);
			return false;
		}
		info.size_axioms[i] = proofs[index];
		info.size_axioms[i]->reference_count++;
		info.size_axioms.length++;
	} for (unsigned int i = 0; i < descendant_count; i++) {
		if (!read(index, in)) {
			free(info);
			return false;
		}
		unsigned int bucket = info.descendants.index_to_insert(index);
		info.descendants.keys[bucket] = index;
		info.descendants.size++;
	} for (unsigned int i = 0; i < element_count; i++) {
		if (!read(info.elements[i], in)) {
			free(info);
			return false;
		}
		info.elements.length++;
	} for (unsigned int i = 0; i < provable_element_count; i++) {
		if (!read(info.provable_elements[i], in)) {
			free(info);
			return false;
		}
		info.provable_elements.length++;
	} for (unsigned int i = 0; i < newly_disjoint_cache_size; i++) {
		if (!read(info.newly_disjoint_cache[i].key, in)
		 || !read(info.newly_disjoint_cache[i].value, in))
		{
			free(info);
			return false;
		}
		info.newly_disjoint_cache.length++;
	}
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Stream>
bool write(
		const set_info<BuiltInConstants, ProofCalculus>& info, Stream& out,
		const hash_map<const typename ProofCalculus::Proof*, unsigned int>& proof_map)
{
	typedef typename ProofCalculus::Proof Proof;

	if (!write(info.arity, out)
	 || !write(info.set_size, out)
	 || !write(info.size_axioms.length, out)
	 || !write(info.descendants.size, out)
	 || !write(info.elements.length, out)
	 || !write(info.provable_elements.length, out)
	 || !write(info.newly_disjoint_cache.length, out))
		return false;
	for (const Proof* proof : info.size_axioms)
		if (!write(proof_map.get(proof), out)) return false;
	for (unsigned int descendant : info.descendants)
		if (!write(descendant, out)) return false;
	for (const tuple_element& element : info.elements)
		if (!write(element, out)) return false;
	for (const tuple& element : info.provable_elements)
		if (!write(element, out)) return false;
	for (const pair<unsigned int, unsigned int>& entry : info.newly_disjoint_cache) {
		if (!write(entry.key, out)
		 || !write(entry.value, out))
			return false;
	}
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int min_set_size, unsigned int max_set_size)
{ }

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool compute_new_set_size(unsigned int set_id,
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int& out, unsigned int min_set_size, unsigned int max_set_size)
{
	out = (max_set_size == UINT_MAX) ? (min_set_size + 200) : ((min_set_size + max_set_size + 1) / 2);
	return true;
}

template<typename Proof>
constexpr bool on_new_size_axiom(Proof* new_size_axiom) {
	return true;
}

template<typename Proof>
inline void on_old_size_axiom(Proof* old_size_axiom) {
	return;
}

struct pair_sorter { };

template<typename K, typename V>
inline bool less_than(const pair<K, V>& first, const pair<K, V>& second, const pair_sorter& sorter) {
	if (first.key < second.key) return true;
	else if (first.key > second.key) return false;
	else if (first.value < second.value) return true;
	else return false;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
struct set_reasoning
{
	typedef typename ProofCalculus::Language Formula;
	typedef typename ProofCalculus::Proof Proof;
	typedef typename Formula::Type FormulaType;

	struct changes {
		array_map<unsigned int, unsigned int> removed_sets;
	};

	extensional_set_graph<ProofCalculus> extensional_graph;
	intensional_set_graph intensional_graph;
	set_info<BuiltInConstants, ProofCalculus>* sets;

	unsigned int capacity;
	unsigned int set_count;

	hash_map<Formula, unsigned int> set_ids;

	hash_multiset<unsigned int> symbols_in_formulas;

	set_reasoning() :
			extensional_graph(1024), intensional_graph(1024),
			capacity(1024), set_count(0), set_ids(2048),
			symbols_in_formulas(256)
	{
		sets = (set_info<BuiltInConstants, ProofCalculus>*) malloc(sizeof(set_info<BuiltInConstants, ProofCalculus>) * capacity);
		if (sets == NULL) exit(EXIT_FAILURE);
		for (unsigned int i = 1; i < capacity; i++)
			sets[i].size_axioms.data = nullptr;

		unsigned int empty_set_id;
		Formula* empty = Formula::new_false();
		if (empty == NULL || !get_set_id_canonicalized(empty, 1, empty_set_id))
			exit(EXIT_FAILURE);
		core::free(*empty); if (empty->reference_count == 0) core::free(empty);
		sets[empty_set_id].change_size(0);
	}

	~set_reasoning() { free_helper(); }

	inline void free_helper() {
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axioms.data != nullptr) {
				extensional_graph.template free_set<false>(i);
				intensional_graph.free_set<false>(i);
				core::free(sets[i]);
			}
		} for (auto entry : set_ids)
			core::free(entry.key);
		core::free(sets);
	}

	static inline void free(set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets) {
		sets.free_helper();
		core::free(sets.extensional_graph);
		core::free(sets.intensional_graph);
		core::free(sets.set_ids);
		core::free(sets.symbols_in_formulas);
	}

	static inline bool clone(
			const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& src,
			set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& dst,
			array_map<const Proof*, Proof*>& proof_map,
			hash_map<const Formula*, Formula*>& formula_map)
	{
		dst.capacity = src.capacity;
		dst.intensional_graph.vertices = (intensional_set_vertex*) malloc(sizeof(intensional_set_vertex) * dst.capacity);
		if (dst.intensional_graph.vertices == nullptr) {
			fprintf(stderr, "set_reasoning.clone ERROR: Insufficient memory for `intensional_set_graph.vertices`.\n");
			return false;
		}
		dst.extensional_graph.vertices = (extensional_set_vertex<ProofCalculus>*) malloc(sizeof(extensional_set_vertex<ProofCalculus>) * dst.capacity);
		if (dst.extensional_graph.vertices == nullptr) {
			fprintf(stderr, "set_reasoning.clone ERROR: Insufficient memory for `extensional_graph.vertices`.\n");
			core::free(dst.intensional_graph.vertices); return false;
		}
		dst.sets = (set_info<BuiltInConstants, ProofCalculus>*) malloc(sizeof(set_info<BuiltInConstants, ProofCalculus>) * dst.capacity);
		if (dst.sets == nullptr) {
			fprintf(stderr, "set_reasoning.clone ERROR: Insufficient memory for `extensional_graph.vertices`.\n");
			core::free(dst.intensional_graph.vertices); core::free(dst.extensional_graph.vertices);
			return false;
		}
		for (unsigned int i = 1; i < dst.capacity; i++)
			dst.sets[i].size_axioms.data = nullptr;
		dst.set_count = src.set_count;
		if (!hash_map_init(dst.set_ids, src.set_ids.table.capacity)) {
			core::free(dst.intensional_graph.vertices);
			core::free(dst.extensional_graph.vertices);
			core::free(dst.sets); return false;
		} else if (!init(dst.symbols_in_formulas, src.symbols_in_formulas.counts.table.capacity)) {
			core::free(dst.intensional_graph.vertices);
			core::free(dst.extensional_graph.vertices);
			core::free(dst.sets); core::free(dst.set_ids);
			return false;
		}

		for (unsigned int i = 1; i < dst.set_count + 1; i++) {
			if (src.sets[i].size_axioms.data != nullptr) {
				if (!extensional_set_vertex<ProofCalculus>::clone_except_parents(src.extensional_graph.vertices[i], dst.extensional_graph.vertices[i], proof_map, formula_map)) {
					core::free(dst);
					return false;
				} else if (!intensional_set_vertex::clone(src.intensional_graph.vertices[i], dst.intensional_graph.vertices[i])) {
					dst.extensional_graph.template free_set<false>(i);
					core::free(dst); return false;
				} else if (!set_info<BuiltInConstants, ProofCalculus>::clone(src.sets[i], dst.sets[i], proof_map, formula_map)) {
					dst.extensional_graph.template free_set<false>(i);
					dst.intensional_graph.template free_set<false>(i);
					core::free(dst); return false;
				}
			}
		} for (unsigned int i = 1; i < dst.set_count + 1; i++) {
			if (src.sets[i].size_axioms.data != nullptr) {
				if (!extensional_set_vertex<ProofCalculus>::clone_only_parents(src.extensional_graph.vertices[i], dst.extensional_graph.vertices[i], proof_map)) {
					core::free(dst);
					return false;
				}
			}
		} for (const auto& entry : src.set_ids) {
			unsigned int bucket = dst.set_ids.table.index_to_insert(entry.key);
			dst.set_ids.values[bucket] = entry.value;
			if (!::clone(entry.key, dst.set_ids.table.keys[bucket], formula_map)) {
				core::free(dst);
				return false;
			}
			dst.set_ids.table.size++;
		} for (const auto& entry : src.symbols_in_formulas.counts) {
			unsigned int bucket = dst.symbols_in_formulas.counts.table.index_to_insert(entry.key);
			dst.symbols_in_formulas.counts.table.keys[bucket] = entry.key;
			dst.symbols_in_formulas.counts.values[bucket] = entry.value;
			dst.symbols_in_formulas.counts.table.size++;
		}
		dst.symbols_in_formulas.sum = src.symbols_in_formulas.sum;
		return true;
	}

	template<bool AncestorsIsEmpty = false>
	inline bool get_ancestors(unsigned int set_id, hash_set<unsigned int>& ancestors) const
	{
		unsigned int index;
		if (AncestorsIsEmpty) {
			index = ancestors.index_to_insert(set_id);
		} else {
			if (!ancestors.check_size()) return false;
			bool contains;
			index = ancestors.index_of(set_id, contains);
			if (contains) return true;
		}

		array<unsigned int> stack(8);
		stack[stack.length++] = set_id;
		ancestors.keys[index] = set_id;
		ancestors.size++;
		while (stack.length != 0) {
			unsigned int current = stack.pop();
			for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (ancestors.contains(parent)) continue;
				if (!ancestors.add(parent) || !stack.add(parent)) {
					return false;
				}
			} for (const auto& entry : extensional_graph.vertices[current].parents) {
				if (ancestors.contains(entry.key)) continue;
				if (!ancestors.add(entry.key) || !stack.add(entry.key)) {
					return false;
				}
			}
		}
		return true;
	}

	template<bool EmptySet>
	bool recompute_provable_elements(unsigned int set_id)
	{
		/* we need to recompute the provable elements of all ancestor nodes; first clear them all */
		array<unsigned int> stack(8);
		hash_set<unsigned int> visited(16);
		stack[stack.length++] = set_id;
		visited.add(set_id);
		while (stack.length > 0) {
			unsigned int current = stack.pop();
			for (tuple& tup : sets[current].provable_elements) core::free(tup);
			sets[current].provable_elements.clear();
			if (!EmptySet || current != set_id) {
				const tuple_element* elements_src = sets[current].elements.data;
				for (unsigned int i = 0; i < sets[current].element_count(); i++) {
					if (!init(sets[current].provable_elements[i], elements_src + (i * sets[current].arity), sets[current].arity))
						return false;
					sets[current].provable_elements.length++;
				}
				if (sets[current].provable_elements.length > 1)
					sort(sets[current].provable_elements, default_sorter());
			}
			for (const auto& entry : extensional_graph.vertices[current].parents) {
				if (visited.contains(entry.key)) continue;
				if (!stack.add(entry.key) || !visited.add(entry.key)) return false;
			} for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (visited.contains(parent)) continue;
				if (!stack.add(parent) || !visited.add(parent)) return false;
			}
		}

		for (unsigned int set : visited)
			if (!stack.add(set)) return false;
		while (stack.length > 0) {
			unsigned int current = stack.pop();
			unsigned int old_provable_element_count = sets[current].provable_elements.length;
			for (unsigned int child : intensional_graph.vertices[current].children) {
				array<tuple> new_provable_elements(max(1, sets[current].provable_elements.length + sets[child].provable_elements.length));
				set_union(new_provable_elements.data, new_provable_elements.length,
						sets[current].provable_elements.data, sets[current].provable_elements.length,
						sets[child].provable_elements.data, sets[child].provable_elements.length);
				swap(new_provable_elements, sets[current].provable_elements);
				for (tuple& tup : new_provable_elements) core::free(tup);
			}
			for (const auto& entry : extensional_graph.vertices[current].children) {
				array<tuple> new_provable_elements(max(1, sets[current].provable_elements.length + sets[entry.key].provable_elements.length));
				set_union(new_provable_elements.data, new_provable_elements.length,
						sets[current].provable_elements.data, sets[current].provable_elements.length,
						sets[entry.key].provable_elements.data, sets[entry.key].provable_elements.length);
				swap(new_provable_elements, sets[current].provable_elements);
				for (tuple& tup : new_provable_elements) core::free(tup);
			}
			if (old_provable_element_count == sets[current].provable_elements.length)
				continue;

			for (const auto& entry : extensional_graph.vertices[current].parents)
				if (!stack.add(entry.key)) return false;
			for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (extensional_graph.vertices[current].parents.contains(parent)) continue;
				if (!stack.add(parent)) return false;
			}
		}
		return true;
	}

	inline void remove_element_at(unsigned int set_id, unsigned int i)
	{
		sets[set_id].remove_element_at(i);
		recompute_provable_elements<false>(set_id);
	}

	bool add_to_provable_elements(unsigned int set_id, const tuple& tup, const hash_set<unsigned int>& ancestors) {
		for (unsigned int ancestor : ancestors) {
			array<tuple>& provable_elements = sets[ancestor].provable_elements;
			unsigned int index = linear_search(provable_elements.data, tup, 0, provable_elements.length);
			if (index == sets[ancestor].provable_elements.length) {
				if (!provable_elements.ensure_capacity(provable_elements.length + 1)
				 || !init(provable_elements[provable_elements.length], tup))
					return false;
				provable_elements.length++;
			} else if (provable_elements[index] != tup) {
				if (!provable_elements.ensure_capacity(provable_elements.length + 1)) return false;
				shift_right(provable_elements.data, provable_elements.length, index);
				if (!init(provable_elements[index], tup))
					return false;
				provable_elements.length++;
			}
		}
		return true;
	}

	inline bool add_to_provable_elements(unsigned int set_id, const tuple& tup) {
		hash_set<unsigned int> ancestors(16);
		if (!get_ancestors<true>(set_id, ancestors))
			return false;
		return add_to_provable_elements(set_id, tup, ancestors);
	}

	inline bool add_element(unsigned int set_id, const tuple& tup, const hash_set<unsigned int>& ancestors) {
		return add_to_provable_elements(set_id, tup, ancestors)
			&& sets[set_id].add_element(tup);
	}

	inline bool add_element(unsigned int set_id, const tuple& tup) {
		hash_set<unsigned int> ancestors(16);
		return get_ancestors<true>(set_id, ancestors)
			&& add_to_provable_elements(set_id, tup, ancestors)
			&& sets[set_id].add_element(tup);
	}

	bool ensure_capacity(unsigned int new_length) {
		if (new_length <= capacity) return true;
		unsigned int new_capacity = capacity;
		expand_capacity(new_capacity, new_length);

		if (!resize(sets, new_capacity) || !extensional_graph.resize(new_capacity) || !intensional_graph.resize(new_capacity))
			return false;
		for (unsigned int i = capacity; i < new_capacity; i++)
			sets[i].size_axioms.data = nullptr; /* we use this field to mark which vertices are free */
		capacity = new_capacity;
		return true;
	}

	unsigned int get_next_free_set() const {
		for (unsigned int i = 1; i < capacity; i++)
			if (sets[i].size_axioms.data == nullptr) return i;
		return capacity;
	}

	bool is_freeable(unsigned int set_id) const {
		if (set_id <= 1) return false;
		for (Proof* size_axiom : sets[set_id].size_axioms) {
			if (size_axiom->reference_count != 1)
				return false;
		}
		return extensional_graph.vertices[set_id].children.size == 0
			&& extensional_graph.vertices[set_id].parents.size == 0;
	}

	template<typename... Args>
	bool is_freeable(unsigned int set_id, array<Proof*>& freeable_axioms, Args&&... visitor) const {
		if (set_id <= 1) return false;
		for (Proof* size_axiom : sets[set_id].size_axioms) {
			if (!freeable_axioms.contains(size_axiom) && size_axiom->reference_count != 1)
				return false;
		}
		return extensional_graph.vertices[set_id].children.size == 0
			&& extensional_graph.vertices[set_id].parents.size == 0;
	}

	template<typename Arg, typename... Args>
	inline bool is_freeable(unsigned int set_id, const Arg& arg, Args&&... visitor) const {
		return is_freeable(set_id, std::forward<Args>(visitor)...);
	}

	template<typename... Args>
	inline void try_free_set(unsigned int set_id, Args&&... visitor) {
		if (is_freeable(set_id, std::forward<Args>(visitor)...)) {
			for (Proof* size_axiom : sets[set_id].size_axioms)
				on_old_size_axiom(size_axiom, std::forward<Args>(visitor)...);
			unsigned int min_set_size; unsigned int max_set_size;
			set_size_bounds(set_id, min_set_size, max_set_size);
			on_free_set(set_id, *this, min_set_size, max_set_size, std::forward<Args>(visitor)...);
			free_set(set_id);
		}
	}

	template<typename... Args>
	inline void try_free_set(Formula* set_formula, unsigned int arity, Args&&... visitor) {
		array_map<unsigned int, unsigned int> variable_map(max(16u, arity));
		for (unsigned int i = 0; i < arity; i++) {
			variable_map.keys[variable_map.size] = i + 1;
			variable_map.values[variable_map.size++] = i + 1;
		}
		Formula* canonicalized = Canonicalizer::canonicalize(*set_formula, variable_map);
		if (canonicalized == nullptr) return;

#if !defined(NDEBUG)
		bool contains;
		unsigned int set_id = set_ids.get(*canonicalized, contains);
		if (!contains)
			fprintf(stderr, "set_reasoning.try_free_set WARNING: The given `set_formula` is not in `set_ids`.\n");
#else
		unsigned int set_id = set_ids.get(*canonicalized);
#endif
		try_free_set(set_id, std::forward<Args>(visitor)...);
		core::free(*canonicalized); if (canonicalized->reference_count == 0) core::free(canonicalized);
	}

	inline bool new_set(unsigned int& set_id)
	{
		set_id = get_next_free_set();
		if (!ensure_capacity(set_id + 1)
		 || !extensional_graph.new_set(set_id))
		{
			return false;
		} else if (!intensional_graph.new_set(set_id)) {
			extensional_graph.template free_set<true>(set_id);
			return false;
		}
		return true;
	}

	template<typename... Args>
	bool new_set(Formula* set_formula, unsigned int arity, unsigned int& set_id, Args&&... visitor)
	{
		if (!new_set(set_id)) return false;

		/* initialize all intensional set relations */
		array<unsigned int> supersets(8), subsets(8);
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axioms.data == nullptr) continue;
			if (sets[i].arity != arity) {
				if (i == 1 && !subsets.add(i)) {
					intensional_graph.template free_set<true>(set_id);
					extensional_graph.template free_set<true>(set_id);
					return false;
				}
				continue;
			}
			if (is_subset(sets[i].set_formula(), set_formula)) {
				if (!subsets.add(i)) {
					intensional_graph.template free_set<true>(set_id);
					extensional_graph.template free_set<true>(set_id);
					return false;
				}
			} else if (is_subset(set_formula, sets[i].set_formula())) {
				if (!supersets.add(i)) {
					intensional_graph.template free_set<true>(set_id);
					extensional_graph.template free_set<true>(set_id);
					return false;
				}
			}
		}

		/* compute the ancestors of `new_set` and their degrees (in the subgraph of the ancestors) */
		array_map<unsigned int, unsigned int> ancestors(8);
		for (unsigned int superset : supersets) {
			if (!::get_ancestors(intensional_graph, superset, ancestors)) {
				intensional_graph.template free_set<true>(set_id);
				extensional_graph.template free_set<true>(set_id);
				return false;
			}
		}
		/* remove ancestors with non-zero out-degree */
		unsigned int next = 0;
		for (unsigned int i = 0; i < ancestors.size; i++) {
			if (ancestors.values[i] == 0) {
				ancestors.keys[next] = ancestors.keys[i];
				ancestors.values[next++] = 0;
			}
		}
		ancestors.size = next;
		/* the vertices with out-degree 0 in `ancestors` are the immediate parents of the new vertex `set_id` */
		for (const auto& entry : ancestors) {
			if (!intensional_graph.add_edge(entry.key, set_id)) {
				intensional_graph.template free_set<true>(set_id);
				extensional_graph.template free_set<true>(set_id);
				return false;
			}
		}

		/* compute the descendants of `new_set` and their degrees (in the subgraph of the descendants) */
		array_map<unsigned int, unsigned int> descendants(8);
		for (unsigned int subset : subsets) {
			if (!get_descendants(intensional_graph, subset, descendants)) {
				intensional_graph.template free_set<true>(set_id);
				extensional_graph.template free_set<true>(set_id);
				return false;
			}
		}
		/* remove descendants with non-zero in-degree */
		next = 0;
		for (unsigned int i = 0; i < descendants.size; i++) {
			if (descendants.values[i] == 0) {
				descendants.keys[next] = descendants.keys[i];
				descendants.values[next++] = 0;
			}
		}
		descendants.size = next;
		/* the vertices with in-degree 0 in `descendants` are the immediate children of the new vertex `set_id` */
		for (const auto& entry : descendants) {
			if (!intensional_graph.add_edge(set_id, entry.key)) {
				intensional_graph.template free_set<true>(set_id);
				extensional_graph.template free_set<true>(set_id);
				return false;
			}
		}

		/* initialize the set_info structure and the set size */
		if (!init(sets[set_id], arity, 1, set_formula, std::forward<Args>(visitor)...)) {
			intensional_graph.template free_set<true>(set_id);
			extensional_graph.template free_set<true>(set_id);
			return false;
		}

		/* remove edges connecting ancestors to descendants */
		for (const auto& immediate_ancestor : ancestors) {
			for (const auto& immediate_descendant : descendants) {
				/* remove the edge from `parent` to `child` if it exists */
				if (intensional_graph.vertices[immediate_ancestor.key].children.contains(immediate_descendant.key))
					intensional_graph.remove_edge(immediate_ancestor.key, immediate_descendant.key);
			}
		}

		/* compute the descendants and provable elements of `set_id` from its immediate children */
		for (unsigned int immediate_descendant : intensional_graph.vertices[set_id].children) {
			if (!sets[set_id].descendants.add_all(sets[immediate_descendant].descendants)) {
				free_set_id(set_id); return false;
			}

			array<tuple> new_provable_elements(max(1, sets[set_id].provable_elements.length + sets[immediate_descendant].provable_elements.length));
			set_union(new_provable_elements.data, new_provable_elements.length,
					sets[set_id].provable_elements.data, sets[set_id].provable_elements.length,
					sets[immediate_descendant].provable_elements.data, sets[immediate_descendant].provable_elements.length);
			swap(new_provable_elements, sets[set_id].provable_elements);
			for (tuple& tup : new_provable_elements) core::free(tup);
		}

		/* update the descendants set of all ancestors of `set_id` */
		hash_set<unsigned int> visited(32);
		array<unsigned int> stack(8);
		stack[stack.length++] = set_id;
		while (stack.length > 0) {
			unsigned int current = stack.pop();
			if (!sets[current].descendants.add(set_id)) {
				free_set_id(set_id); return false;
			}

			for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (visited.contains(parent)) continue;
				if (!visited.add(parent) || !stack.add(parent)) {
					free_set_id(set_id); return false;
				}
			} for (const auto& entry : extensional_graph.vertices[current].parents) {
				if (visited.contains(entry.key)) continue;
				if (!visited.add(entry.key) || !stack.add(entry.key)) {
					free_set_id(set_id); return false;
				}
			}
		}

		/* precompute entries for `newly_disjoint_cache` */
		array<unsigned int> strictly_partial_subsets(8);
		array<unsigned int> strictly_partial_supersets(8);
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (i == set_id || sets[i].size_axioms.data == nullptr || sets[i].arity != arity)
				continue;
			if (is_strictly_partial_subset(sets[i].set_formula(), set_formula)) {
				if (!strictly_partial_subsets.add(i)) {
					free_set_id(set_id); return false;
				}
			} if (is_strictly_partial_subset(set_formula, sets[i].set_formula())) {
				if (!strictly_partial_supersets.add(i)) {
					free_set_id(set_id); return false;
				}
			}
		}
		for (unsigned int i = 0; i < strictly_partial_subsets.length; i++) {
			unsigned int first = strictly_partial_subsets[i];
			for (unsigned int j = i + 1; j < strictly_partial_subsets.length; j++) {
				unsigned int second = strictly_partial_subsets[j];
				if (are_newly_disjoint(sets[first].set_formula(), sets[second].set_formula(), set_id)) {
					if (!sets[set_id].newly_disjoint_cache.add(make_pair(first, second))) {
						free_set_id(set_id); return false;
					}
				}
			}
		}
		/* ensure that the `newly_disjoint_cache` of the parents of the new
		   node don't contain any entries that are also in the
		   `newly_disjoint_cache` of this node */
		for (unsigned int parent : intensional_graph.vertices[set_id].parents) {
			for (unsigned int i = 0; i < sets[parent].newly_disjoint_cache.length; i++) {
				if (sets[set_id].newly_disjoint_cache.contains(sets[parent].newly_disjoint_cache[i]))
					sets[parent].newly_disjoint_cache.remove(i--);
			}
		} for (const auto& entry : extensional_graph.vertices[set_id].parents) {
			for (unsigned int i = 0; i < sets[entry.key].newly_disjoint_cache.length; i++) {
				if (sets[set_id].newly_disjoint_cache.contains(sets[entry.key].newly_disjoint_cache[i]))
					sets[entry.key].newly_disjoint_cache.remove(i--);
			}
		}
		for (unsigned int partial_superset : strictly_partial_supersets) {
			strictly_partial_subsets.clear();
			for (unsigned int i = 1; i < set_count + 1; i++) {
				if (i == partial_superset || i == set_id || sets[i].size_axioms.data == nullptr || sets[i].arity != arity)
					continue;
				if (is_strictly_partial_subset(sets[i].set_formula(), sets[partial_superset].set_formula())) {
					if (!strictly_partial_subsets.add(i)) {
						free_set_id(set_id); return false;
					}
				}
			}
			for (unsigned int first : strictly_partial_subsets) {
				if (are_newly_disjoint(sets[first].set_formula(), set_formula, partial_superset)) {
					if (!sets[partial_superset].newly_disjoint_cache.add(first < set_id ? make_pair(first, set_id) : make_pair(set_id, first))) {
						free_set_id(set_id); return false;
					}
				}
			}
		}

		/* compute the upper bound and lower bound on the size of this new set */
		unsigned int min_set_size, max_set_size, initial_set_size;
		if (!set_size_bounds(set_id, min_set_size, max_set_size)
		 || !compute_new_set_size(set_id, *this, initial_set_size, min_set_size, max_set_size, std::forward<Args>(visitor)...))
		{
			free_set_id(set_id); return false;
		}
		sets[set_id].change_size(initial_set_size);

		if (set_id == set_count + 1)
			set_count++;
		return true;
	}

	inline bool set_size_bounds(unsigned int set_id,
			unsigned int& min_set_size, unsigned int& max_set_size) const
	{
		if (!is_unfixed(set_id)) {
			min_set_size = sets[set_id].set_size;
			max_set_size = min_set_size;
			return true;
		}

		max_set_size = UINT_MAX;
		for (unsigned int parent : intensional_graph.vertices[set_id].parents) {
			if (sets[parent].set_size == 0) {
				max_set_size = 0;
				break;
			}
		}
		return get_size_lower_bound(set_id, min_set_size)
			&& (max_set_size == 0 || get_size_upper_bound(set_id, max_set_size));
	}

	bool free_set_id(unsigned int set_id)
	{
#if !defined(NDEBUG)
		if (extensional_graph.vertices[set_id].parents.size > 0
		 || extensional_graph.vertices[set_id].children.size > 0)
			fprintf(stderr, "set_reasoning.free_set WARNING: This set has extensional edges.\n");
#endif

		for (unsigned int parent : intensional_graph.vertices[set_id].parents) {
			for (const pair<unsigned int, unsigned int>& newly_disjoint_pair : sets[set_id].newly_disjoint_cache) {
				if (are_newly_disjoint(sets[newly_disjoint_pair.key].set_formula(), sets[newly_disjoint_pair.value].set_formula(), parent))
					sets[parent].newly_disjoint_cache.add(newly_disjoint_pair);
			}
		} for (const auto& entry : extensional_graph.vertices[set_id].parents) {
			for (const pair<unsigned int, unsigned int>& newly_disjoint_pair : sets[set_id].newly_disjoint_cache) {
				if (are_newly_disjoint(sets[newly_disjoint_pair.key].set_formula(), sets[newly_disjoint_pair.value].set_formula(), entry.key))
					sets[entry.key].newly_disjoint_cache.add(newly_disjoint_pair);
			}
		}
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (i == set_id || sets[i].size_axioms.data == nullptr) continue;
			for (unsigned int j = 0; j < sets[i].newly_disjoint_cache.length; j++) {
				if (sets[i].newly_disjoint_cache[j].key == set_id || sets[i].newly_disjoint_cache[j].value == set_id)
					sets[i].newly_disjoint_cache.remove(j--);
			}
		}

		recompute_provable_elements<true>(set_id);

		hash_set<unsigned int> visited(32);
		array<unsigned int> stack(8);
		stack[stack.length++] = set_id;
		while (stack.length > 0) {
			unsigned int current = stack.pop();
			sets[current].descendants.remove(set_id);

			for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (visited.contains(parent)) continue;
				if (!visited.add(parent) || !stack.add(parent))
					return false;
			} for (const auto& entry : extensional_graph.vertices[current].parents) {
				if (visited.contains(entry.key)) continue;
				if (!visited.add(entry.key) || !stack.add(entry.key))
					return false;
			}
		}

		const array<unsigned int>& parents_src = intensional_graph.vertices[set_id].parents;
		const array<unsigned int>& children_src = intensional_graph.vertices[set_id].children;
		array<unsigned int> parents(max((size_t) 1, parents_src.length));
		array<unsigned int> children(max((size_t) 1, children_src.length));
		for (unsigned int i = 0; i < parents_src.length; i++)
			parents[i] = parents_src[i];
		parents.length = parents_src.length;
		for (unsigned int i = 0; i < children_src.length; i++)
			children[i] = children_src[i];
		children.length = children_src.length;

		intensional_graph.free_set<true>(set_id);
		extensional_graph.template free_set<true>(set_id);

		for (unsigned int parent : parents) {
			array_map<unsigned int, unsigned int> descendants(8);
			if (!get_descendants(intensional_graph, parent, descendants))
				return false;
			for (unsigned int child : children) {
				/* first check if there is already a path from `parent` to `child` */
				if (descendants.contains(child)) continue;
				if (!intensional_graph.add_edge(parent, child)) return false;
			}
		}

		core::free(sets[set_id]);
		if (set_id == set_count + 2) set_count--;
		return true;
	}

	inline bool free_set(unsigned int set_id)
	{
		array_multiset<unsigned int> symbols(16);
		Formula* formula = sets[set_id].set_formula();
		if (!get_constants(*formula, symbols)) return false;
		symbols_in_formulas.subtract<true>(symbols);

		bool contains;
		unsigned int bucket = set_ids.table.index_of(*formula, contains);
		core::free(set_ids.table.keys[bucket]);
		set_ids.remove_at(bucket);
		return free_set_id(set_id);
	}

	template<typename... Args>
	inline bool get_set_id_canonicalized(Formula* formula, unsigned int arity, unsigned int& set_id, bool& is_new, Args&&... visitor) {
		bool contains; unsigned int bucket;
		if (!set_ids.check_size()) return false;
		set_id = set_ids.get(*formula, contains, bucket);
		if (!contains) {
			array_multiset<unsigned int> symbols(16);
			if (!get_constants(*formula, symbols)) return false;
			symbols_in_formulas.add(symbols);

			if (!init(set_ids.table.keys[bucket], *formula)) {
				return false;
			} else if (!new_set(formula, arity, set_id, std::forward<Args>(visitor)...)) {
				core::free(set_ids.table.keys[bucket]);
				core::set_empty(set_ids.table.keys[bucket]);
				return false;
			}
			set_ids.values[bucket] = set_id;
			set_ids.table.size++;
			is_new = true;
		} else {
			is_new = false;
		}
		return true;
	}

	template<typename... Args>
	inline bool get_set_id(Formula* formula, unsigned int arity, unsigned int& set_id, bool& is_new, Args&&... visitor) {
		array_map<unsigned int, unsigned int> variable_map(max(16u, arity));
		for (unsigned int i = 0; i < arity; i++) {
			variable_map.keys[variable_map.size] = i + 1;
			variable_map.values[variable_map.size++] = i + 1;
		}
		Formula* canonicalized = Canonicalizer::canonicalize(*formula, variable_map);
		if (canonicalized == nullptr)
			return false;
		bool result = get_set_id_canonicalized(canonicalized, arity, set_id, is_new, std::forward<Args>(visitor)...);
		core::free(*canonicalized); if (canonicalized->reference_count == 0) core::free(canonicalized);
		return result;
	}

	template<typename... Args>
	inline bool get_set_id_canonicalized(Formula* formula, unsigned int arity, unsigned int& set_id, Args&&... visitor) {
		bool is_new;
		return get_set_id_canonicalized(formula, arity, set_id, is_new, std::forward<Args>(visitor)...);
	}

	template<typename... Args>
	inline bool get_set_id(Formula* formula, unsigned int arity, unsigned int& set_id, Args&&... visitor) {
		bool is_new;
		return get_set_id(formula, arity, set_id, is_new, std::forward<Args>(visitor)...);
	}

	inline unsigned int get_existing_set_id(Formula* formula, unsigned int arity) const
	{
		array_map<unsigned int, unsigned int> variable_map(max(16u, arity));
		for (unsigned int i = 0; i < arity; i++) {
			variable_map.keys[variable_map.size] = i + 1;
			variable_map.values[variable_map.size++] = i + 1;
		}
		Formula* canonicalized_formula = Canonicalizer::canonicalize(*formula, variable_map);
		if (canonicalized_formula == nullptr) return UINT_MAX;

#if !defined(NDEBUG)
		bool contains;
		unsigned int set_id = set_ids.get(*canonicalized_formula, contains);
		if (!contains)
			fprintf(stderr, "set_reasoning.get_existing_set_id WARNING: No such set for given formula.\n");
#else
		unsigned int set_id = set_ids.get(*canonicalized_formula);
#endif
		core::free(*canonicalized_formula); if (canonicalized_formula->reference_count == 0) core::free(canonicalized_formula);
		return set_id;
	}

	inline bool get_strongly_connected_component(
			unsigned int set, hash_set<unsigned int>& component)
	{
		array<unsigned int> stack(8);
		component.add(set);
		stack[stack.length++] = set;
		while (stack.length > 0) {
			unsigned int current = stack.pop();

			for (const auto& entry : extensional_graph.vertices[current].children) {
				if (!sets[entry.key].descendants.contains(set) || component.contains(entry.key)) continue;
				if (!component.add(entry.key) || !stack.add(entry.key)) return false;
			} for (unsigned int child : intensional_graph.vertices[current].children) {
				if (!sets[child].descendants.contains(set) || component.contains(child)) continue;
				if (!component.add(child) || !stack.add(child)) return false;
			}
		}
		return true;
	}

	inline bool get_strongly_connected_component(
			unsigned int set, hash_set<unsigned int>& component,
			const pair<unsigned int, unsigned int>& ignore_edge)
	{
		array<unsigned int> stack(8);
		array<unsigned int> dfs(8);
		stack[stack.length++] = set;
		dfs[dfs.length++] = set;
		while (stack.length > 0) {
			unsigned int current = stack.pop();

			for (const auto& entry : extensional_graph.vertices[current].parents) {
				if ((current == ignore_edge.key && entry.key == ignore_edge.value) || dfs.contains(entry.key)) continue;
				if (!dfs.add(entry.key) || !stack.add(entry.key)) return false;
			} for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (dfs.contains(parent)) continue;
				if (!dfs.add(parent) || !stack.add(parent)) return false;
			}
		}

		array_map<unsigned int, unsigned int> roots(2 * dfs.length);
		for (unsigned int i = 0; i < dfs.length; i++) {
			roots.keys[i] = dfs[i];
			roots.values[i] = 0;
		}
		roots.size = dfs.length;
		for (unsigned int current : dfs) {
			stack[stack.length++] = current;
			while (stack.length > 0) {
				unsigned int inner = stack.pop();
				bool contains;
				unsigned int& root = roots.get(inner, contains);
				if (!contains || root != 0) continue;
				root = current;

				for (const auto& entry : extensional_graph.vertices[current].children) {
					if (entry.key == ignore_edge.key && inner == ignore_edge.value) continue;
					if (!stack.add(entry.key)) return false;
				} for (unsigned int child : intensional_graph.vertices[current].children) {
					if (!stack.add(child)) return false;
				}
			}
		}

		unsigned int root = roots.get(set);
		for (unsigned int i = 0; i < roots.size; i++) {
			if (roots.values[i] == root)
				component.add(roots.keys[i]);
		}
		return true;
	}

	struct subgraph_formula_view {
		set_info<BuiltInConstants, ProofCalculus>* sets;
		unsigned int* indices;
		unsigned int length;

		subgraph_formula_view(
				set_info<BuiltInConstants, ProofCalculus>* sets,
				const hash_set<unsigned int>& vertices) : sets(sets), length(vertices.size)
		{
			indices = (unsigned int*) malloc(sizeof(unsigned int) * vertices.size);
			unsigned int index = 0;
			for (unsigned int vertex : vertices)
				indices[index++] = vertex;
		}

		~subgraph_formula_view() { core::free(indices); }

		Formula* operator[] (unsigned int index) const {
			return sets[indices[index]].set_formula();
		}

		inline unsigned int size() const {
			return length;
		}
	};

	bool uncontract_component(unsigned int contracted_set,
			const hash_set<unsigned int>& connected_component,
			const array<pair<unsigned int, pair<unsigned int, unsigned int>>>& old_disjoint_cache_items)
	{
		/* restore the original `newly_disjoint_cache` */
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axioms.data == nullptr) continue;
			for (unsigned int j = 0; j < sets[i].newly_disjoint_cache.length; j++) {
				const pair<unsigned int, unsigned int>& entry = sets[i].newly_disjoint_cache[j];
				if (entry.key == contracted_set || entry.value == contracted_set)
					sets[i].newly_disjoint_cache.remove(j--);
			}
		}
		for (const pair<unsigned int, pair<unsigned int, unsigned int>>& old_entry : old_disjoint_cache_items)
			sets[old_entry.key].newly_disjoint_cache[sets[old_entry.key].newly_disjoint_cache.length++] = old_entry.value;

		hash_set<unsigned int> visited(16);
		array<unsigned int> stack(8);
		stack[stack.length++] = contracted_set;
		while (stack.length > 0) {
			unsigned int ancestor = stack.pop();
			for (unsigned int member : connected_component)
				sets[ancestor].descendants.add(member);
			sets[ancestor].descendants.remove(contracted_set);

			for (unsigned int parent : intensional_graph.vertices[ancestor].parents) {
				if (visited.contains(parent)) continue;
				if (!visited.add(parent) || !stack.add(parent)) return false;
			} for (const auto& entry : extensional_graph.vertices[ancestor].parents) {
				if (visited.contains(entry.key)) continue;
				if (!visited.add(entry.key) || !stack.add(entry.key)) return false;
			}
		}

		for (unsigned int member : connected_component) {
			for (const auto& entry : extensional_graph.vertices[member].parents) {
				if (!connected_component.contains(entry.key)) {
					array_map<unsigned int, array<Proof*>>& children = extensional_graph.vertices[entry.key].children;
					children.remove(contracted_set);
					array<Proof*>& src = extensional_graph.vertices[member].parents.get(entry.key);
					children.keys[children.size] = member;
					if (!array_init(children.values[children.size], 1 << (core::log2(src.length) + 1)))
						return false;
					for (Proof* proof : src)
						children.values[children.size][children.values[children.size].length++] = proof;
					children.size++;
				}
			} for (unsigned int parent : intensional_graph.vertices[member].parents) {
				if (!connected_component.contains(parent)) {
					unsigned int index = intensional_graph.vertices[parent].children.index_of(contracted_set);
					if (index < intensional_graph.vertices[parent].children.length)
						intensional_graph.vertices[parent].children.remove(index);
					if (!intensional_graph.vertices[parent].children.contains(member))
						intensional_graph.vertices[parent].children.add(member);
				}
			} for (const auto& entry : extensional_graph.vertices[member].children) {
				if (!connected_component.contains(entry.key)) {
					array_map<unsigned int, array<Proof*>>& parents = extensional_graph.vertices[entry.key].parents;
					parents.remove(contracted_set);
					array<Proof*>& src = extensional_graph.vertices[member].children.get(entry.key);
					parents.keys[parents.size] = member;
					if (!array_init(parents.values[parents.size], 1 << (core::log2(src.length) + 1)))
						return false;
					for (Proof* proof : src)
						parents.values[parents.size][parents.values[parents.size].length++] = proof;
					parents.size++;
				}
			} for (unsigned int child : intensional_graph.vertices[member].children) {
				if (!connected_component.contains(child)) {
					unsigned int index = intensional_graph.vertices[child].parents.index_of(contracted_set);
					if (index < intensional_graph.vertices[child].parents.length)
						intensional_graph.vertices[child].parents.remove(index);
					if (!intensional_graph.vertices[child].parents.contains(member))
						intensional_graph.vertices[child].parents.add(member);
				}
			}
		}

		intensional_graph.vertices[contracted_set].parents.clear();
		intensional_graph.vertices[contracted_set].children.clear();
		extensional_graph.vertices[contracted_set].parents.clear();
		extensional_graph.vertices[contracted_set].children.clear();

		intensional_graph.template free_set<true>(contracted_set);
		extensional_graph.template free_set<true>(contracted_set);
		core::free(sets[contracted_set]); if (contracted_set == set_count + 2) set_count--;
		return true;
	}

	bool contract_component(
			unsigned int set, unsigned int& contracted_set,
			const hash_set<unsigned int>& connected_component, bool& is_fixed,
			array<pair<unsigned int, pair<unsigned int, unsigned int>>>& old_disjoint_cache_items)
	{
		is_fixed = false;

		Formula* conjunction = Formula::new_and(subgraph_formula_view(sets, connected_component));
		if (conjunction == NULL) return false;
		for (unsigned int member : connected_component)
			sets[member].set_formula()->reference_count++;

		if (!new_set(contracted_set)) {
			core::free(*conjunction); if (conjunction->reference_count == 0) core::free(conjunction);
			return false;
		} else if (!init(sets[contracted_set], sets[set].arity, sets[set].set_size, conjunction)) {
			core::free(*conjunction); if (conjunction->reference_count == 0) core::free(conjunction);
			intensional_graph.template free_set<true>(contracted_set);
			extensional_graph.template free_set<true>(contracted_set);
			return false;
		}
		core::free(*conjunction); if (conjunction->reference_count == 0) core::free(conjunction);
		if (contracted_set == set_count + 1) set_count++;

		if (!sets[contracted_set].descendants.add_all(sets[set].descendants)) {
			intensional_graph.template free_set<true>(contracted_set);
			extensional_graph.template free_set<true>(contracted_set);
			core::free(sets[contracted_set]); if (contracted_set == set_count + 2) set_count--;
			return false;
		}
		for (const tuple& tup : sets[set].provable_elements) {
			if (!sets[contracted_set].provable_elements.add(tup)) {
				intensional_graph.template free_set<true>(contracted_set);
				extensional_graph.template free_set<true>(contracted_set);
				core::free(sets[contracted_set]); if (contracted_set == set_count + 2) set_count--;
				return false;
			}
		}

		for (unsigned int member : connected_component) {
			if (!is_unfixed(member))
				is_fixed = true;
			for (const auto& entry : extensional_graph.vertices[member].parents) {
				if (!connected_component.contains(entry.key)) {
					array_map<unsigned int, array<Proof*>>& children = extensional_graph.vertices[entry.key].children;
					array_map<unsigned int, array<Proof*>>& parents = extensional_graph.vertices[contracted_set].parents;
					unsigned int index = children.index_of(member);
					core::free(children.values[index]);
					children.remove_at(index);
					children.keys[children.size++] = contracted_set;
					if (!parents.ensure_capacity(parents.size + 1)) {
						uncontract_component(contracted_set, connected_component, old_disjoint_cache_items);
						return false;
					}
					parents.keys[parents.size++] = entry.key;
				}
			} for (unsigned int parent : intensional_graph.vertices[member].parents) {
				if (!connected_component.contains(parent)) {
					unsigned int index = intensional_graph.vertices[parent].children.index_of(member);
					intensional_graph.vertices[parent].children.remove(index);
					if (!intensional_graph.vertices[parent].children.contains(contracted_set))
						intensional_graph.vertices[parent].children.add(contracted_set);
					if (!intensional_graph.vertices[contracted_set].parents.contains(parent)
					 && !intensional_graph.vertices[contracted_set].parents.add(parent)) {
						uncontract_component(contracted_set, connected_component, old_disjoint_cache_items);
						return false;
					}
				}
			} for (const auto& entry : extensional_graph.vertices[member].children) {
				if (!connected_component.contains(entry.key)) {
					array_map<unsigned int, array<Proof*>>& parents = extensional_graph.vertices[entry.key].parents;
					array_map<unsigned int, array<Proof*>>& children = extensional_graph.vertices[contracted_set].children;
					unsigned int index = parents.index_of(member);
					core::free(parents.values[index]);
					parents.remove_at(index);
					parents.keys[parents.size++] = contracted_set;
					if (!children.ensure_capacity(children.size + 1)) {
						uncontract_component(contracted_set, connected_component, old_disjoint_cache_items);
						return false;
					}
					children.keys[children.size++] = entry.key;
				}
			} for (unsigned int child : intensional_graph.vertices[member].children) {
				if (!connected_component.contains(child)) {
					unsigned int index = intensional_graph.vertices[child].parents.index_of(member);
					intensional_graph.vertices[child].parents.remove(index);
					if (!intensional_graph.vertices[child].parents.contains(contracted_set))
						intensional_graph.vertices[child].parents.add(contracted_set);
					if (!intensional_graph.vertices[contracted_set].children.contains(child)
					 && !intensional_graph.vertices[contracted_set].children.add(child)) {
						uncontract_component(contracted_set, connected_component, old_disjoint_cache_items);
						return false;
					}
				}
			}
		}

		hash_set<unsigned int> visited(16);
		array<unsigned int> stack(8);
		stack[stack.length++] = contracted_set;
		while (stack.length > 0) {
			unsigned int ancestor = stack.pop();
			for (unsigned int member : connected_component)
				sets[ancestor].descendants.remove(member);
			sets[ancestor].descendants.add(contracted_set);

			for (unsigned int parent : intensional_graph.vertices[ancestor].parents) {
				if (visited.contains(parent)) continue;
				if (!visited.add(parent) || !stack.add(parent)) {
					uncontract_component(contracted_set, connected_component, old_disjoint_cache_items);
					return false;
				}
			} for (const auto& entry : extensional_graph.vertices[ancestor].parents) {
				if (visited.contains(entry.key)) continue;
				if (!visited.add(entry.key) || !stack.add(entry.key)) {
					uncontract_component(contracted_set, connected_component, old_disjoint_cache_items);
					return false;
				}
			}
		}

		/* add the elements of the sets in the component to the new set */
		for (unsigned int set : connected_component) {
			if (!sets[contracted_set].elements.ensure_capacity(sets[contracted_set].elements.length + sets[set].elements.length)) {
				uncontract_component(contracted_set, connected_component, old_disjoint_cache_items);
				return false;
			}
			for (const tuple_element& src_element : sets[set].elements) {
				tuple_element& new_element = sets[contracted_set].elements[sets[contracted_set].elements.length];
				if (!init(new_element, src_element)) {
					uncontract_component(contracted_set, connected_component, old_disjoint_cache_items);
					return false;
				}
				sets[contracted_set].elements.length++;
			}
		}

		/* compute `newly_disjoint_cache` for the new set */
		for (unsigned int set : connected_component) {
			if (!sets[contracted_set].newly_disjoint_cache.append(sets[set].newly_disjoint_cache.data, sets[set].newly_disjoint_cache.length)) {
				uncontract_component(contracted_set, connected_component, old_disjoint_cache_items);
				return false;
			}
		}
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axioms.data == nullptr) continue;
			if (!old_disjoint_cache_items.ensure_capacity(old_disjoint_cache_items.length + sets[i].newly_disjoint_cache.length)) {
				uncontract_component(contracted_set, connected_component, old_disjoint_cache_items);
				return false;
			}
			for (unsigned int j = 0; j < sets[i].newly_disjoint_cache.length; j++) {
				const pair<unsigned int, unsigned int>& entry = sets[i].newly_disjoint_cache[j];
				if (connected_component.contains(entry.key)) {
					old_disjoint_cache_items[old_disjoint_cache_items.length++] = make_pair(i, entry);
					sets[i].newly_disjoint_cache.remove(j--);
					if (!connected_component.contains(entry.value))
						sets[i].newly_disjoint_cache[sets[i].newly_disjoint_cache.length++] = make_pair(contracted_set, entry.value);
				} else if (connected_component.contains(entry.value)) {
					old_disjoint_cache_items[old_disjoint_cache_items.length++] = make_pair(i, entry);
					sets[i].newly_disjoint_cache.remove(j--);
					sets[i].newly_disjoint_cache[sets[i].newly_disjoint_cache.length++] = make_pair(entry.key, contracted_set);
				}
			}
		}
		return true;
	}

	bool increase_set_size(
			unsigned int set, unsigned int requested_size,
			array<unsigned int>& stack, bool& graph_changed)
	{
		/* first check if `set` is part of a strongly-connected component */
		bool is_fixed;
		hash_set<unsigned int> connected_component(16);
		array<pair<unsigned int, pair<unsigned int, unsigned int>>> old_disjoint_cache_items(32);
		if (!get_strongly_connected_component(set, connected_component)) {
			return false;
		} else if (connected_component.size > 1) {
			/* perform "surgery" on the graph by creating a new "virtual"
			   vertex that represents the strongly-connected component with all
			   of its internal edges contracted */
			if (!contract_component(set, set, connected_component, is_fixed, old_disjoint_cache_items))
				return false;
		} else {
			is_fixed = !is_unfixed(set);
		}

		if (!stack.ensure_capacity(stack.length + 1)) {
			if (connected_component.size > 1)
				uncontract_component(set, connected_component, old_disjoint_cache_items);
			return false;
		} if (stack.contains(set)) {
			/* we have previously requested to change the size of `set` */
			graph_changed = false; return true;
		}
		stack[stack.length++] = set;

		unsigned int* clique = nullptr; unsigned int clique_count, ancestor_of_clique, upper_bound = 0;
		if (is_fixed) {
			upper_bound = sets[set].set_size;
		} else if (!get_size_upper_bound(set, upper_bound, clique, clique_count, ancestor_of_clique)) {
			if (connected_component.size > 1)
				uncontract_component(set, connected_component, old_disjoint_cache_items);
			return false;
		}
		if (requested_size <= upper_bound) {
			/* the bound is already satisfied */
			core::free(clique); stack.length--;
			graph_changed = true;
			if (connected_component.size > 1) {
				for (unsigned int member : connected_component)
					sets[member].change_size(requested_size);
				return uncontract_component(set, connected_component, old_disjoint_cache_items);
			} else {
				sets[set].change_size(requested_size);
				return true;
			}
		} else if (sets[set].set_size < upper_bound) {
			core::free(clique); stack.length--;
			graph_changed = true;
			if (connected_component.size > 1) {
				for (unsigned int member : connected_component)
					sets[member].change_size(upper_bound);
				return uncontract_component(set, connected_component, old_disjoint_cache_items);
			} else {
				sets[set].change_size(upper_bound);
				return true;
			}
		}

		if (!is_fixed && clique != nullptr) {
			bool child_graph_changed = false;
			unsigned int clique_upper_bound = sets[ancestor_of_clique].set_size;
			for (unsigned int i = 0; i < clique_count; i++)
				clique_upper_bound -= sets[clique[i]].set_size;
			if (clique_upper_bound < requested_size && !increase_set_size(ancestor_of_clique, sets[ancestor_of_clique].set_size + (requested_size - upper_bound), stack, child_graph_changed))
			{
				if (connected_component.size > 1)
					uncontract_component(set, connected_component, old_disjoint_cache_items);
				core::free(clique); return false;
			} else if (child_graph_changed) {
				/* the graph has been changed */
				core::free(clique); stack.length--;
				graph_changed = true;
				if (connected_component.size > 1)
					return uncontract_component(set, connected_component, old_disjoint_cache_items);
				else return true;
			}

			for (unsigned int i = 0; i < clique_count; i++) {
				unsigned int requested_child_size;
				if (sets[clique[i]].set_size > requested_size - upper_bound)
					requested_child_size = sets[clique[i]].set_size - (requested_size - upper_bound);
				else requested_child_size = 0;
				if (!decrease_set_size(clique[i], requested_child_size, stack, child_graph_changed)) {
					if (connected_component.size > 1)
						uncontract_component(set, connected_component, old_disjoint_cache_items);
					core::free(clique); return false;
				} else if (child_graph_changed) {
					/* the graph has been changed */
					core::free(clique); stack.length--;
					graph_changed = true;
					if (connected_component.size > 1)
						uncontract_component(set, connected_component, old_disjoint_cache_items);
					else return true;
				}
			}
		}

		/* we were unable to change any bounds */
		core::free(clique); stack.length--;
		graph_changed = false;
		if (connected_component.size > 1)
			return uncontract_component(set, connected_component, old_disjoint_cache_items);
		else return true;
	}

	bool decrease_set_size(
			unsigned int set, unsigned int requested_size,
			array<unsigned int>& stack, bool& graph_changed)
	{
		/* first check if `set` is part of a strongly-connected component */
		bool is_fixed;
		hash_set<unsigned int> connected_component(16);
		array<pair<unsigned int, pair<unsigned int, unsigned int>>> old_disjoint_cache_items(32);
		if (!get_strongly_connected_component(set, connected_component)) {
			return false;
		} else if (connected_component.size > 1) {
			/* perform "surgery" on the graph by creating a new "virtual"
			   vertex that represents the strongly-connected component with all
			   of its internal edges contracted */
			if (!contract_component(set, set, connected_component, is_fixed, old_disjoint_cache_items))
				return false;
		} else {
			is_fixed = !is_unfixed(set);
		}

		if (!stack.ensure_capacity(stack.length + 1)) {
			if (connected_component.size > 1)
				uncontract_component(set, connected_component, old_disjoint_cache_items);
			return false;
		} if (stack.contains(set)) {
			/* we have previously requested to change the size of `set` */
			graph_changed = false; return true;
		}
		stack[stack.length++] = set;

		unsigned int* clique = NULL; unsigned int clique_count, lower_bound, ancestor_of_clique;
		if (is_fixed) {
			lower_bound = sets[set].set_size;
		} else if (!get_size_lower_bound(set, lower_bound, clique, clique_count, ancestor_of_clique)) {
			if (connected_component.size > 1)
				uncontract_component(set, connected_component, old_disjoint_cache_items);
			return false;
		}
		if (requested_size >= lower_bound) {
			/* the bound is already satisfied */
			core::free(clique); stack.length--;
			graph_changed = true;
			if (connected_component.size > 1) {
				for (unsigned int member : connected_component)
					sets[member].change_size(requested_size);
				return uncontract_component(set, connected_component, old_disjoint_cache_items);
			} else {
				sets[set].change_size(requested_size);
				return true;
			}
		} else if (sets[set].set_size > lower_bound) {
			core::free(clique); stack.length--;
			graph_changed = true;
			if (connected_component.size > 1) {
				for (unsigned int member : connected_component)
					sets[member].change_size(lower_bound);
				return uncontract_component(set, connected_component, old_disjoint_cache_items);
			} else {
				sets[set].change_size(lower_bound);
				return true;
			}
		}

		if (!is_fixed) {
			bool child_graph_changed;
			if (ancestor_of_clique == 0) {
				/* the lower bound comes from a subset clique of `set` */
				for (unsigned int i = 0; i < clique_count; i++) {
					unsigned int requested_child_size;
					if (sets[clique[i]].set_size > lower_bound - requested_size)
						requested_child_size = sets[clique[i]].set_size - (lower_bound - requested_size);
					else requested_child_size = 0;
					if (!decrease_set_size(clique[i], requested_child_size, stack, child_graph_changed)) {
						if (connected_component.size > 1)
							uncontract_component(set, connected_component, old_disjoint_cache_items);
						core::free(clique); return false;
					} else if (child_graph_changed) {
						/* the graph has been changed */
						core::free(clique); stack.length--;
						graph_changed = true;
						if (connected_component.size > 1)
							return uncontract_component(set, connected_component, old_disjoint_cache_items);
						else return true;
					}
				}
			} else {
				/* the lower bound comes from an ancestor of `set` and the `requested_size` is 0 */
				unsigned int requested_ancestor_size = 0;
				for (unsigned int i = 0; i < clique_count; i++)
					requested_ancestor_size += sets[clique[i]].set_size;
				if (!increase_set_size(ancestor_of_clique, requested_ancestor_size, stack, child_graph_changed)) {
					if (connected_component.size > 1)
						uncontract_component(set, connected_component, old_disjoint_cache_items);
					core::free(clique); return false;
				} else if (child_graph_changed) {
					/* the graph has been changed */
					core::free(clique); stack.length--;
					graph_changed = true;
					if (connected_component.size > 1)
						return uncontract_component(set, connected_component, old_disjoint_cache_items);
					else return true;
				}

				for (unsigned int i = 0; i < clique_count; i++) {
					unsigned int requested_child_size;
					if (sets[clique[i]].set_size > requested_ancestor_size - sets[ancestor_of_clique].set_size)
						requested_child_size = sets[clique[i]].set_size - (requested_ancestor_size - sets[ancestor_of_clique].set_size);
					else requested_child_size = 0;
					if (!decrease_set_size(clique[i], requested_child_size, stack, child_graph_changed)) {
						if (connected_component.size > 1)
							uncontract_component(set, connected_component, old_disjoint_cache_items);
						core::free(clique); return false;
					} else if (child_graph_changed) {
						/* the graph has been changed */
						core::free(clique); stack.length--;
						graph_changed = true;
						if (connected_component.size > 1)
							return uncontract_component(set, connected_component, old_disjoint_cache_items);
						else return true;
					}
				}
			}
		}

		/* we were unable to change any bounds */
		core::free(clique); stack.length--;
		graph_changed = false;
		if (connected_component.size > 1)
			return uncontract_component(set, connected_component, old_disjoint_cache_items);
		else return true;
	}

	template<typename... Args>
	bool add_existing_subset_axiom(
		Proof* axiom, Formula* antecedent, Formula* consequent,
		unsigned int antecedent_set, unsigned int consequent_set)
	{
#if !defined(NDEBUG)
		if (consequent_set == antecedent_set)
			fprintf(stderr, "set_reasoning.add_subset_axiom WARNING: `consequent` and `antecedent` are the same set.\n");
#endif

		if (!extensional_graph.add_edge(consequent_set, antecedent_set, axiom))
			return false;

		if (!update_descendants(antecedent_set, consequent_set)) {
			remove_subset_relation<true>(antecedent_set, consequent_set, antecedent, consequent);
			return false;
		}

		/* update the provable elements of `consequent_set` and all its ancestors */
		array<unsigned int> ancestors_stack(8);
		hash_set<unsigned int> visited(16);
		ancestors_stack[ancestors_stack.length++] = consequent_set;
		visited.add(consequent_set);
		visited.add(antecedent_set);
		while (ancestors_stack.length > 0) {
			unsigned int current = ancestors_stack.pop();
			array<tuple> new_provable_elements(max(1, sets[antecedent_set].provable_elements.length + sets[current].provable_elements.length));
			set_union(new_provable_elements.data, new_provable_elements.length,
					sets[antecedent_set].provable_elements.data, sets[antecedent_set].provable_elements.length,
					sets[current].provable_elements.data, sets[current].provable_elements.length);
			swap(new_provable_elements, sets[current].provable_elements);
			for (tuple& tup : new_provable_elements) core::free(tup);

			for (const auto& entry : extensional_graph.vertices[current].parents) {
				if (visited.contains(entry.key)) continue;
				if (!ancestors_stack.add(entry.key) || !visited.add(entry.key)) {
					extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
					return false;
				}
			} for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (visited.contains(parent)) continue;
				if (!ancestors_stack.add(parent) || !visited.add(parent)) {
					extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
					return false;
				}
			}
		}
		return true;
	}

	/* NOTE: this function does not check for consistency */
	template<typename... Args>
	bool add_subset_axiom(Proof* axiom, Args&&... visitor)
	{
		Formula* operand = axiom->formula->quantifier.operand;
		unsigned int arity = 1;
		while (operand->type == FormulaType::FOR_ALL) {
			operand = operand->quantifier.operand;
			arity++;
		}
		Formula* antecedent = operand->binary.left;
		Formula* consequent = operand->binary.right;

		if (!set_ids.check_size(set_ids.table.size + 2)) return false;

		unsigned int antecedent_set, consequent_set;
		bool is_antecedent_new, is_consequent_new;
		if (!get_set_id(antecedent, arity, antecedent_set, is_antecedent_new, std::forward<Args>(visitor)...)) {
			return false;
		} else if (!get_set_id(consequent, arity, consequent_set, is_consequent_new, std::forward<Args>(visitor)...)) {
			if (is_freeable(antecedent_set)) free_set(antecedent_set);
			return false;
		}

		if (!add_existing_subset_axiom(axiom, antecedent, consequent, antecedent_set, consequent_set)) {
			/* if either the antecedent or consequent sets have no references, free them */
			if (is_freeable(consequent_set)) {
				for (Proof* size_axiom : sets[consequent_set].size_axioms)
					on_old_size_axiom(size_axiom, std::forward<Args>(visitor)...);
				on_free_set(consequent_set, *this, 0, 0);
				free_set(consequent_set);
			} if (is_freeable(antecedent_set)) {
				for (Proof* size_axiom : sets[antecedent_set].size_axioms)
					on_old_size_axiom(size_axiom, std::forward<Args>(visitor)...);
				on_free_set(antecedent_set, *this, 0, 0);
				free_set(antecedent_set);
			}
			return false;
		}
		return true;
	}

	template<bool ResolveInconsistencies, typename... Args>
	Proof* get_existing_subset_axiom(Formula* antecedent, Formula* consequent, unsigned int arity, Args&&... visitor) const
	{
		array_map<unsigned int, unsigned int> variable_map(max(16u, arity));
		for (unsigned int i = 0; i < arity; i++) {
			variable_map.keys[variable_map.size] = i + 1;
			variable_map.values[variable_map.size++] = i + 1;
		}
		Formula* canonicalized_antecedent = Canonicalizer::canonicalize(*antecedent, variable_map);
		if (canonicalized_antecedent == nullptr) return nullptr;
		variable_map.size = arity;
		Formula* canonicalized_consequent = Canonicalizer::canonicalize(*consequent, variable_map);
		if (canonicalized_consequent == nullptr) {
			core::free(*canonicalized_antecedent); if (canonicalized_antecedent->reference_count == 0) core::free(canonicalized_antecedent);
			return nullptr;
		}

#if !defined(NDEBUG)
		bool contains;
		unsigned int antecedent_set = set_ids.get(*canonicalized_antecedent, contains);
		if (!contains)
			fprintf(stderr, "set_reasoning.get_existing_subset_axiom WARNING: No such set for given antecedent.\n");
		unsigned int consequent_set = set_ids.get(*canonicalized_consequent, contains);
		if (!contains)
			fprintf(stderr, "set_reasoning.get_existing_subset_axiom WARNING: No such set for given consequent.\n");
#else
		unsigned int antecedent_set = set_ids.get(*canonicalized_antecedent);
		unsigned int consequent_set = set_ids.get(*canonicalized_consequent);
#endif
		core::free(*canonicalized_antecedent); if (canonicalized_antecedent->reference_count == 0) core::free(canonicalized_antecedent);
		core::free(*canonicalized_consequent); if (canonicalized_consequent->reference_count == 0) core::free(canonicalized_consequent);

#if !defined(NDEBUG)
		if (consequent_set == antecedent_set)
			fprintf(stderr, "set_reasoning.get_existing_subset_axiom WARNING: `consequent` and `antecedent` are the same set.\n");
#endif

		return extensional_graph.get_existing_edge(consequent_set, antecedent_set, consequent, antecedent);
	}

	inline bool update_descendants(unsigned int antecedent_set, unsigned int consequent_set)
	{
		/* update the descendants of `consequent_set` and all its ancestors */
		array<pair<unsigned int, unsigned int>> stack(8);
		stack[stack.length++] = {consequent_set, antecedent_set};
		while (stack.length > 0) {
			pair<unsigned int, unsigned int> current = stack.pop();
			unsigned int old_size = sets[current.key].descendants.size;
			if (!sets[current.key].descendants.add_all(sets[current.value].descendants))
				return false;
			if (old_size == sets[current.key].descendants.size) continue;

			for (unsigned int parent : intensional_graph.vertices[current.key].parents) {
				if (!stack.add({parent, current.key}))
					return false;
			} for (const auto& entry : extensional_graph.vertices[current.key].parents) {
				if (!stack.add({entry.key, current.key}))
					return false;
			}
		}
		return true;
	}

	inline unsigned int count_union_of_provable_elements(unsigned int first, unsigned int second) const {
		unsigned int count = sets[first].provable_elements.length + sets[second].provable_elements.length;
		auto union_both = [&count](const tuple& tup, unsigned int i, unsigned int j) { count--; };
		auto union_one = [](const tuple& tup, unsigned int i, unsigned int j) { };
		set_union(union_both, union_one, union_one,
				sets[first].provable_elements.data, sets[first].provable_elements.length,
				sets[second].provable_elements.data, sets[second].provable_elements.length);
		return count;
	}

	template<bool ResolveInconsistencies, bool FreeSets, typename... Args>
	Proof* get_subset_axiom(
			Formula* antecedent, Formula* consequent, unsigned int arity,
			unsigned int& antecedent_set, unsigned int& consequent_set,
			bool& is_antecedent_new, bool& is_consequent_new, Args&&... visitor)
	{
		if (!set_ids.check_size(set_ids.table.size + 2)) return nullptr;

		if (!get_set_id(antecedent, arity, antecedent_set, is_antecedent_new, std::forward<Args>(visitor)...)) {
			return nullptr;
		} else if (!get_set_id(consequent, arity, consequent_set, is_consequent_new, std::forward<Args>(visitor)...)) {
			if (FreeSets) try_free_set(antecedent_set);
			return nullptr;
		}

#if !defined(NDEBUG)
		if (consequent_set == antecedent_set)
			fprintf(stderr, "set_reasoning.get_subset_axiom WARNING: `consequent` and `antecedent` are the same set.\n");
#endif

		bool new_edge;
		Proof* axiom = extensional_graph.get_edge(consequent_set, antecedent_set, consequent, antecedent, arity, new_edge);
		if (axiom == nullptr) {
			/* if either the antecedent or consequent sets have no references, free them */
			if (FreeSets) {
				try_free_set(consequent_set);
				try_free_set(antecedent_set);
			}
			return nullptr;
		}

		if (!new_edge) return axiom;

		/* check that the new edge does not create any inconsistencies */
		array<unsigned int> ancestors_stack(8);
		hash_set<unsigned int> visited(16);
		ancestors_stack[ancestors_stack.length++] = consequent_set;
		visited.add(consequent_set);
		visited.add(antecedent_set);
		while (ancestors_stack.length > 0) {
			unsigned int current = ancestors_stack.pop();
			unsigned int provable_element_count = count_union_of_provable_elements(antecedent_set, current);
			while (provable_element_count > sets[current].set_size) {
				if (ResolveInconsistencies) {
					bool graph_changed;
					array<unsigned int> stack(8);
					if (!increase_set_size(current, provable_element_count, stack, graph_changed) || !graph_changed) {
						extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
						if (FreeSets) {
							try_free_set(consequent_set);
							try_free_set(antecedent_set);
						}
						return nullptr;
					}
				} else {
					extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
					if (FreeSets) {
						try_free_set(consequent_set);
						try_free_set(antecedent_set);
					}
					return nullptr;
				}
			}

			for (const auto& entry : extensional_graph.vertices[current].parents) {
				if (visited.contains(entry.key)) continue;
				if (!ancestors_stack.add(entry.key) || !visited.add(entry.key)) {
					extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
					if (FreeSets) {
						try_free_set(consequent_set);
						try_free_set(antecedent_set);
					}
					return nullptr;
				}
			} for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (visited.contains(parent)) continue;
				if (!ancestors_stack.add(parent) || !visited.add(parent)) {
					extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
					if (FreeSets) {
						try_free_set(consequent_set);
						try_free_set(antecedent_set);
					}
					return nullptr;
				}
			}
		}
		for (unsigned int descendant : sets[antecedent_set].descendants) {
			if (sets[descendant].set_size > 0 && are_disjoint(consequent_set, descendant)) {
				if (ResolveInconsistencies) {
					while (sets[descendant].set_size > 0) {
						bool graph_changed;
						array<unsigned int> stack(8);
						if (!decrease_set_size(descendant, 0, stack, graph_changed)) {
							extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
							if (FreeSets) {
								try_free_set(consequent_set);
								try_free_set(antecedent_set);
							}
							return nullptr;
						} else if (graph_changed) continue;

						/* we were unable to change the graph, so it cannot be made consistent (as far as we can tell) */
						extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
						if (FreeSets) {
							try_free_set(consequent_set);
							try_free_set(antecedent_set);
						}
						return nullptr;
					}
				} else {
					extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
					if (FreeSets) {
						try_free_set(consequent_set);
						try_free_set(antecedent_set);
					}
					return nullptr;
				}
			}
		}
		while (true) {
			unsigned int* clique = nullptr; unsigned int clique_count; unsigned int ancestor_of_clique;
			if (sets[antecedent_set].set_size == 0) {
				break;
			} else if (sets[consequent_set].set_size == 0) {
				bool graph_changed;
				array<unsigned int> stack(8);
				if (!decrease_set_size(antecedent_set, 0, stack, graph_changed)) {
					extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
					if (FreeSets) {
						try_free_set(consequent_set);
						try_free_set(antecedent_set);
					}
					return nullptr;
				} else if (graph_changed) continue;

				/* we were unable to change the graph, so it cannot be made consistent (as far as we can tell) */
				extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
				if (FreeSets) {
					try_free_set(consequent_set);
					try_free_set(antecedent_set);
				}
				return nullptr;
			} else if (!find_largest_disjoint_clique_with_set(*this, antecedent_set, consequent_set, clique, clique_count, ancestor_of_clique, 1))
			{
				extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
				if (FreeSets) {
					try_free_set(consequent_set);
					try_free_set(antecedent_set);
				}
			}

			if (clique == nullptr) break;

			if (ResolveInconsistencies) {
				/* try increasing the size of the ancestor */
				unsigned int requested_size = 0;
				for (unsigned int i = 0; i < clique_count; i++)
					requested_size += sets[clique[i]].set_size;

				bool graph_changed;
				array<unsigned int> stack(8);
				if (!increase_set_size(ancestor_of_clique, requested_size, stack, graph_changed)) {
					extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
					if (FreeSets) {
						try_free_set(consequent_set);
						try_free_set(antecedent_set);
					}
					core::free(clique); return nullptr;
				} else if (graph_changed) { core::free(clique); continue; }

				for (unsigned int i = 0; i < clique_count; i++) {
					unsigned int new_requested_size;
					if (sets[ancestor_of_clique].set_size > requested_size - sets[clique[i]].set_size)
						new_requested_size = sets[ancestor_of_clique].set_size - (requested_size - sets[clique[i]].set_size);
					else new_requested_size = 0;

					if (!decrease_set_size(clique[i], new_requested_size, stack, graph_changed)) {
						extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
						if (FreeSets) {
							try_free_set(consequent_set);
							try_free_set(antecedent_set);
						}
						core::free(clique); return nullptr;
					} else if (graph_changed) break;
				}

				if (graph_changed) { core::free(clique); continue; }

				/* we were unable to change the graph, so it cannot be made consistent (as far as we can tell) */
				extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
				if (FreeSets) {
					try_free_set(consequent_set);
					try_free_set(antecedent_set);
				}
				core::free(clique); return nullptr;
			} else {
				extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
				if (FreeSets) {
					try_free_set(consequent_set);
					try_free_set(antecedent_set);
				}
				core::free(clique); return nullptr;
			}
		}

		/* update the descendants of `consequent_set` and all its ancestors */
		if (!update_descendants(antecedent_set, consequent_set)) {
			remove_subset_relation<true>(antecedent_set, consequent_set, antecedent, consequent);
			return nullptr;
		}

		/* update the provable elements of `consequent_set` and all its ancestors */
		ancestors_stack[ancestors_stack.length++] = consequent_set;
		visited.clear();
		visited.add(consequent_set);
		visited.add(antecedent_set);
		while (ancestors_stack.length > 0) {
			unsigned int current = ancestors_stack.pop();
			array<tuple> new_provable_elements(max(1, sets[antecedent_set].provable_elements.length + sets[current].provable_elements.length));
			set_union(new_provable_elements.data, new_provable_elements.length,
					sets[antecedent_set].provable_elements.data, sets[antecedent_set].provable_elements.length,
					sets[current].provable_elements.data, sets[current].provable_elements.length);
			swap(new_provable_elements, sets[current].provable_elements);
			for (tuple& tup : new_provable_elements) core::free(tup);

			for (const auto& entry : extensional_graph.vertices[current].parents) {
				if (visited.contains(entry.key)) continue;
				if (!ancestors_stack.add(entry.key) || !visited.add(entry.key)) {
					extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
					if (FreeSets) {
						try_free_set(consequent_set);
						try_free_set(antecedent_set);
					}
					return nullptr;
				}
			} for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (visited.contains(parent)) continue;
				if (!ancestors_stack.add(parent) || !visited.add(parent)) {
					extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);
					if (FreeSets) {
						try_free_set(consequent_set);
						try_free_set(antecedent_set);
					}
					return nullptr;
				}
			}
		}
		return axiom;
	}

	template<bool FreeSets, typename... Args>
	bool remove_subset_relation(
			unsigned int antecedent_set, unsigned int consequent_set,
			Formula* antecedent, Formula* consequent, Args&&... visitor)
	{
#if !defined(NDEBUG)
		if (consequent_set == antecedent_set)
			fprintf(stderr, "set_reasoning.remove_subset_relation WARNING: `consequent` and `antecedent` are the same set.\n");
#endif

		extensional_graph.remove_edge(consequent_set, antecedent_set, consequent, antecedent);

		/* we need to recompute the descendants and provable elements of all ancestor nodes; first clear them all */
		array<unsigned int> stack(8);
		hash_set<unsigned int> visited(16);
		stack[stack.length++] = consequent_set;
		visited.add(consequent_set);
		while (stack.length > 0) {
			unsigned int current = stack.pop();
			sets[current].descendants.clear();
			sets[current].descendants.add(current);
			for (tuple& tup : sets[current].provable_elements) core::free(tup);
			sets[current].provable_elements.clear();
			const tuple_element* elements_src = sets[current].elements.data;
			for (unsigned int i = 0; i < sets[current].element_count(); i++) {
				if (!init(sets[current].provable_elements[i], elements_src + (i * sets[current].arity), sets[current].arity))
					return false;
				sets[current].provable_elements.length++;
			}
			if (sets[current].provable_elements.length > 1)
				sort(sets[current].provable_elements, default_sorter());
			for (const auto& entry : extensional_graph.vertices[current].parents) {
				if (visited.contains(entry.key)) continue;
				if (!stack.add(entry.key) || !visited.add(entry.key)) return false;
			} for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (visited.contains(parent)) continue;
				if (!stack.add(parent) || !visited.add(parent)) return false;
			}
		}

		stack[stack.length++] = consequent_set;
		while (stack.length > 0) {
			unsigned int current = stack.pop();
			unsigned int old_descendant_count = sets[current].descendants.size;
			unsigned int old_provable_element_count = sets[current].provable_elements.length;
			for (unsigned int child : intensional_graph.vertices[current].children) {
				if (!sets[current].descendants.add_all(sets[child].descendants)) return false;
				array<tuple> new_provable_elements(max(1, sets[current].provable_elements.length + sets[child].provable_elements.length));
				set_union(new_provable_elements.data, new_provable_elements.length,
						sets[current].provable_elements.data, sets[current].provable_elements.length,
						sets[child].provable_elements.data, sets[child].provable_elements.length);
				swap(new_provable_elements, sets[current].provable_elements);
				for (tuple& tup : new_provable_elements) core::free(tup);
			}
			for (const auto& entry : extensional_graph.vertices[current].children) {
				if (!sets[current].descendants.add_all(sets[entry.key].descendants)) return false;
				array<tuple> new_provable_elements(max(1, sets[current].provable_elements.length + sets[entry.key].provable_elements.length));
				set_union(new_provable_elements.data, new_provable_elements.length,
						sets[current].provable_elements.data, sets[current].provable_elements.length,
						sets[entry.key].provable_elements.data, sets[entry.key].provable_elements.length);
				swap(new_provable_elements, sets[current].provable_elements);
				for (tuple& tup : new_provable_elements) core::free(tup);
			}
			if (old_descendant_count == sets[current].descendants.size
			 && old_provable_element_count == sets[current].provable_elements.length)
				continue;

			for (const auto& entry : extensional_graph.vertices[current].parents)
				if (!stack.add(entry.key)) return false;
			for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (extensional_graph.vertices[current].parents.contains(parent)) continue;
				if (!stack.add(parent)) return false;
			}
		}

		/* if either the antecedent or consequent sets have no references, free them */
		if (FreeSets) {
			try_free_set(consequent_set, std::forward<Args>(visitor)...);
			try_free_set(antecedent_set, std::forward<Args>(visitor)...);
		}
		return true;
	}

	template<bool FreeSets, typename... Args>
	bool free_subset_axiom(Formula* antecedent, Formula* consequent, unsigned int arity, Args&&... visitor)
	{
		array_map<unsigned int, unsigned int> variable_map(max(16u, arity));
		for (unsigned int i = 0; i < arity; i++) {
			variable_map.keys[variable_map.size] = i + 1;
			variable_map.values[variable_map.size++] = i + 1;
		}
		Formula* canonicalized_antecedent = Canonicalizer::canonicalize(*antecedent, variable_map);
		if (canonicalized_antecedent == nullptr) return false;
		variable_map.size = arity;
		Formula* canonicalized_consequent = Canonicalizer::canonicalize(*consequent, variable_map);
		if (canonicalized_consequent == nullptr) {
			core::free(*canonicalized_antecedent); if (canonicalized_antecedent->reference_count == 0) core::free(canonicalized_antecedent);
			return false;
		}

#if !defined(NDEBUG)
		bool contains;
		unsigned int antecedent_set = set_ids.get(*canonicalized_antecedent, contains);
		if (!contains)
			fprintf(stderr, "set_reasoning.try_free_subset_axiom WARNING: No such set for given antecedent.\n");
		unsigned int consequent_set = set_ids.get(*canonicalized_consequent, contains);
		if (!contains)
			fprintf(stderr, "set_reasoning.try_free_subset_axiom WARNING: No such set for given consequent.\n");
#else
		unsigned int antecedent_set = set_ids.get(*canonicalized_antecedent);
		unsigned int consequent_set = set_ids.get(*canonicalized_consequent);
#endif
		core::free(*canonicalized_antecedent); if (canonicalized_antecedent->reference_count == 0) core::free(canonicalized_antecedent);
		core::free(*canonicalized_consequent); if (canonicalized_consequent->reference_count == 0) core::free(canonicalized_consequent);

		return remove_subset_relation<FreeSets>(antecedent_set, consequent_set, antecedent, consequent, std::forward<Args>(visitor)...);
	}

	template<typename... Args>
	inline bool free_subset_axiom(Proof* subset_axiom, Args&&... visitor)
	{
		Formula* operand = subset_axiom->formula->quantifier.operand;
		unsigned int arity = 1;
		while (operand->type == FormulaType::FOR_ALL) {
			operand = operand->quantifier.operand;
			arity++;
		}
		Formula* antecedent = operand->binary.left;
		Formula* consequent = operand->binary.right;
		return free_subset_axiom<true>(antecedent, consequent, arity, std::forward<Args>(visitor)...);
	}

	template<bool ResolveInconsistencies>
	bool set_size_axiom(unsigned int set_id, unsigned int new_size)
	{
		bool is_fixed = !is_unfixed(set_id);
		if (new_size > sets[set_id].set_size) {
			unsigned int upper_bound = 0;
			if (is_fixed) {
				upper_bound = sets[set_id].set_size;
			} else if (!get_size_upper_bound(set_id, upper_bound))
				return false;
			if (new_size <= upper_bound) {
				sets[set_id].change_size(new_size);
				return true;
			}

			if (!ResolveInconsistencies || is_fixed) return false;

			while (true) {
				/* compute the upper bound on the size of this set; if the new size
				   violates this bound, change the sizes of other sets to increase the bound */
				array<unsigned int> stack(8); bool graph_changed;
				if (!increase_set_size(set_id, new_size, stack, graph_changed))
					return false;
				if (sets[set_id].set_size == new_size)
					break;
				if (!graph_changed) return false;
			}
		} else if (new_size < sets[set_id].set_size) {
			unsigned int lower_bound = 0;
			if (is_fixed) {
				lower_bound = sets[set_id].set_size;
			} else if (!get_size_lower_bound(set_id, lower_bound))
				return false;
			if (new_size >= lower_bound) {
				sets[set_id].change_size(new_size);
				return true;
			}

			if (!ResolveInconsistencies || is_fixed) return false;

			while (true) {
				/* compute the lower bound on the size of this set; if the new size
				   violates this bound, change the sizes of other sets to increase the bound */
				array<unsigned int> stack(8); bool graph_changed;
				if (!decrease_set_size(set_id, new_size, stack, graph_changed))
					return false;
				if (sets[set_id].set_size == new_size)
					break;
				if (!graph_changed) return false;
			}
		}
		return true;
	}

	template<bool ResolveInconsistencies, typename... Args>
	inline Proof* get_size_axiom(
			Formula* formula, unsigned int arity, unsigned int new_size,
			unsigned int& set_id, bool& is_set_new, Args&&... visitor)
	{
		if (!get_set_id(formula, arity, set_id, is_set_new, std::forward<Args>(visitor)...)
		 || !set_size_axiom<ResolveInconsistencies>(set_id, new_size))
			return nullptr;
		return sets[set_id].get_size_axiom(formula, std::forward<Args>(visitor)...);
	}

	bool get_size_lower_bound(unsigned int set_id, unsigned int& lower_bound,
			unsigned int*& clique, unsigned int& clique_count,
			unsigned int& ancestor_of_clique) const
	{
#if !defined(NDEBUG)
		if (!is_unfixed(set_id))
			fprintf(stderr, "get_size_lower_bound WARNING: Set with ID %u is fixed. This function should not be called on fixed sets.\n", set_id);
#endif

		unsigned int intensional_lower_bound = 0;
		Formula* formula = sets[set_id].set_formula();
		if (formula->type == FormulaType::EQUALS && formula->binary.right->type == FormulaType::CONSTANT) {
			lower_bound = 1;
			return true;
		} else if (formula->type == FormulaType::OR) {
			for (unsigned int i = 0; i < formula->array.length; i++) {
				Formula* operand = formula->array.operands[i];
				if (operand->type == FormulaType::EQUALS && operand->binary.right->type == FormulaType::CONSTANT) {
					intensional_lower_bound = max(intensional_lower_bound, 1);
					continue;
				}

				unsigned int subset_id; bool contains;
				subset_id = set_ids.get(*operand, contains);
				if (contains)
					intensional_lower_bound = max(intensional_lower_bound, sets[subset_id].set_size);
			}
		}

		/* first compute the number of elements that provably belong to this set */
		lower_bound = sets[set_id].provable_elements.length;

		/* next compute the maximal clique of subsets of this set
		   (maximal in the sense of the size of their union) */
		if (!find_largest_disjoint_subset_clique(*this, set_id, clique, clique_count, lower_bound))
			return false;

		if (clique != NULL) {
			lower_bound = 0;
			for (unsigned int i = 0; i < clique_count; i++)
				lower_bound += sets[clique[i]].set_size;
		}

		if (lower_bound == 0) {
			/* we need to be more careful since this set being empty causes
			   some ancestor sets to become disjoint, which could cause set
			   size bounds of some ancestors to be violated */

			/* determine which pairs of sets are "newly disjoint" if this `set_id` becomes empty */
			unsigned int old_size = sets[set_id].set_size;
			sets[set_id].set_size = 0;
			for (const pair<unsigned int, unsigned int>& newly_disjoint_pair : sets[set_id].newly_disjoint_cache) {
				/* check if any ancestor of `i` and `j` has a set size smaller than its lower bound */
				if (!find_largest_disjoint_clique_with_edge(*this, newly_disjoint_pair.key, newly_disjoint_pair.value, clique, clique_count, ancestor_of_clique, INT_MIN, 0)) {
					sets[set_id].set_size = old_size;
					return false;
				}
				if (clique != NULL) {
					unsigned int clique_size = 0;
					for (unsigned int i = 0; i < clique_count; i++)
						clique_size += sets[clique[i]].set_size;
					if (clique_size > sets[ancestor_of_clique].set_size) {
						lower_bound = 1;
						break;
					} else {
						core::free(clique);
						clique = NULL;
						clique_count = 0;
						ancestor_of_clique = 0;
					}
				}
			}
			sets[set_id].set_size = old_size;
		} else {
			ancestor_of_clique = 0;
		}
		lower_bound = max(intensional_lower_bound, lower_bound);
		return true;
	}

	inline bool get_size_lower_bound(unsigned int set_id, unsigned int& lower_bound) const
	{
		unsigned int* clique = NULL; unsigned int clique_count, ancestor_of_clique;
		if (!get_size_lower_bound(set_id, lower_bound, clique, clique_count, ancestor_of_clique))
			return false;
		if (clique != NULL) core::free(clique);
		return true;
	}

	bool get_size_upper_bound(unsigned int set_id, unsigned int& upper_bound,
			unsigned int*& clique, unsigned int& clique_count, unsigned int& ancestor_of_clique) const
	{
#if !defined(NDEBUG)
		if (!is_unfixed(set_id))
			fprintf(stderr, "get_size_upper_bound WARNING: Set with ID %u is fixed. This function should not be called on fixed sets.\n", set_id);
#endif

		unsigned int intensional_upper_bound = UINT_MAX;
		Formula* formula = sets[set_id].set_formula();
		if (formula->type == FormulaType::EQUALS && formula->binary.right->type == FormulaType::CONSTANT) {
			upper_bound = 1;
			return true;
		} else if (formula->type == FormulaType::OR) {
			intensional_upper_bound = 0;
			for (unsigned int i = 0; i < formula->array.length; i++) {
				Formula* operand = formula->array.operands[i];
				if (operand->type == FormulaType::EQUALS && operand->binary.right->type == FormulaType::CONSTANT) {
					intensional_upper_bound += 1;
					continue;
				}

				unsigned int subset_id; bool contains;
				subset_id = set_ids.get(*operand, contains);
				if (!contains) {
					intensional_upper_bound = UINT_MAX;
					break;
				} else {
					intensional_upper_bound += sets[subset_id].set_size;
				}
			}
		} else if (formula->type == FormulaType::AND) {
			for (unsigned int i = 0; i < formula->array.length; i++) {
				Formula* operand = formula->array.operands[i];
				if (operand->type == FormulaType::EQUALS && operand->binary.right->type == FormulaType::CONSTANT) {
					upper_bound = 1;
					return true;
				}

				unsigned int subset_id; bool contains;
				subset_id = set_ids.get(*operand, contains);
				if (contains)
					intensional_upper_bound = min(sets[subset_id].set_size, intensional_upper_bound);
			}
		}

		if (!find_largest_disjoint_clique_with_set(*this, set_id, clique, clique_count, ancestor_of_clique, INT_MIN))
			return false;
		if (clique == NULL) {
			upper_bound = intensional_upper_bound; /* the upper bound is infinite */
			return true;
		}

		unsigned int index = index_of(set_id, clique, clique_count);
		clique[index] = clique[--clique_count];

		upper_bound = sets[ancestor_of_clique].set_size;
		for (unsigned int i = 0; i < clique_count; i++)
			upper_bound -= sets[clique[i]].set_size;
		upper_bound = min(upper_bound, intensional_upper_bound);
		return true;
	}

	inline bool get_size_upper_bound(unsigned int set_id, unsigned int& upper_bound) const
	{
		unsigned int* clique = NULL; unsigned int clique_count; unsigned int ancestor_of_clique;
		if (!get_size_upper_bound(set_id, upper_bound, clique, clique_count, ancestor_of_clique))
			return false;
		if (clique != NULL) core::free(clique);
		return true;
	}

	inline bool are_disjoint(Formula* first, Formula* second) const
	{
		Formula* intersection = intersect(first, second);
		if (intersection == nullptr)
			return true;

		array<unsigned int> stack(8);
		hash_set<unsigned int> visited(16);
		stack[stack.length++] = 1;
		visited.add(1);
		while (stack.length > 0) {
			unsigned int current = stack.pop();
			if (is_subset(intersection, sets[current].set_formula())) {
				core::free(*intersection); if (intersection->reference_count == 0) core::free(intersection);
				return true;
			}

			for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (visited.contains(parent) || sets[parent].set_size != 0) continue;
				if (!visited.add(parent) || !stack.add(parent)) exit(EXIT_FAILURE);
			} for (const auto& entry : extensional_graph.vertices[current].parents) {
				if (visited.contains(entry.key) || sets[entry.key].set_size != 0) continue;
				if (!visited.add(entry.key) || !stack.add(entry.key)) exit(EXIT_FAILURE);
			}
		}
		core::free(*intersection); if (intersection->reference_count == 0) core::free(intersection);
		return false;
	}

	inline bool are_disjoint(unsigned int first_set, unsigned int second_set) const
	{
		if (sets[first_set].arity != sets[second_set].arity)
			return false;
		return are_disjoint(sets[first_set].set_formula(), sets[second_set].set_formula());
	}

	/**
	 * Returns `true` when `first` and `second` are disjoint iff the set with
	 * ID `set_id` is empty (and all other sets stay the same).
	 */
	inline bool are_newly_disjoint(Formula* first, Formula* second, unsigned int set_id) const
	{
		Formula* intersection = intersect(first, second);
		if (intersection == nullptr) return false;

		/* first check if the intersection belongs to `set_id` */
		if (!is_subset(intersection, sets[set_id].set_formula())) {
			core::free(*intersection); if (intersection->reference_count == 0) core::free(intersection);
			return false;
		}

		/* then check if the intersection belongs to any child */
		for (const auto& entry : extensional_graph.vertices[set_id].children) {
			if (is_subset(intersection, sets[entry.key].set_formula())) {
				core::free(*intersection); if (intersection->reference_count == 0) core::free(intersection);
				return false;
			}
		} for (unsigned int child : intensional_graph.vertices[set_id].children) {
			if (is_subset(intersection, sets[child].set_formula())) {
				core::free(*intersection); if (intersection->reference_count == 0) core::free(intersection);
				return false;
			}
		}

		core::free(*intersection); if (intersection->reference_count == 0) core::free(intersection);
		return true;
	}

	inline bool is_unfixed(Proof* size_axiom) const {
		return (size_axiom->children.length == 0);
	}

	inline bool is_unfixed(Proof* size_axiom, const array<Proof*>& observations) const {
		return (size_axiom->children.length == 0
			&& !observations.contains(size_axiom));
	}

	inline bool is_unfixed(unsigned int set_id) const {
		for (Proof* size_axiom : sets[set_id].size_axioms)
			if (!is_unfixed(size_axiom)) return false;
		return true;
	}

	inline bool is_unfixed(unsigned int set_id, const array<Proof*>& observations) const {
		for (Proof* size_axiom : sets[set_id].size_axioms)
			if (!is_unfixed(size_axiom, observations)) return false;
		return true;
	}

	bool get_unfixed_sets(array<unsigned int>& unfixed_set_ids, const array<Proof*>& observations) {
		for (unsigned int i = 2; i < set_count + 1; i++) {
			if (sets[i].size_axioms.data == nullptr || !is_unfixed(i, observations))
				continue;
			if (!unfixed_set_ids.add(i)) return false;
		}
		return true;
	}

	bool compute_descendants(unsigned int set, hash_set<unsigned int>& descendants) const {
		array<unsigned int> stack(8);
		stack[stack.length++] = set;
		if (!descendants.add(set)) return false;
		while (stack.length > 0) {
			unsigned int current = stack.pop();

			for (unsigned int child : intensional_graph.vertices[current].children) {
				if (descendants.contains(child)) continue;
				if (!descendants.add(child) || !stack.add(child)) return false;
			} for (const auto& entry : extensional_graph.vertices[current].children) {
				if (descendants.contains(entry.key)) continue;
				if (!descendants.add(entry.key) || !stack.add(entry.key)) return false;
			}
		}
		return true;
	}

	bool are_descendants_valid() const {
		bool success = true;
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axioms.data != nullptr) {
				hash_set<unsigned int> descendants(16);
				if (!compute_descendants(i, descendants)) return false;
				if (!descendants.equals(sets[i].descendants)) {
					print("set_reasoning.are_descendants_valid WARNING: Actual `descendants` doesn't match expected `descendants` for set with ID ", stderr);
					print(i, stderr); print(".\n", stderr);
					print("  Computed: ", stderr); print(descendants, stderr); print('\n', stderr);
					print("  Expected: ", stderr); print(sets[i].descendants, stderr); print('\n', stderr);
					success = false;
				}
			}
		}
		return success;
	}

	template<bool PrintProvableElements = false, typename Stream, typename... Printer>
	bool print_axioms(Stream& out, Printer&&... printer) const {
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axioms.data == nullptr) continue;
			if (!print("Set ID ", out) || !print(i, out) || !print(" has set size axioms [", out))
				return false;
			bool first = true;
			for (Proof* size_axiom : sets[i].size_axioms) {
				if (!first && !print(", ", out)) return false;
				if (!print(*size_axiom->formula, out, std::forward<Printer>(printer)...))
					return false;
				first = false;
			}
			if (!print("]\n", out)) return false;
			for (const auto& entry : extensional_graph.vertices[i].children) {
				for (const Proof* proof : entry.value) {
					if (!print("  ", out) || !print(*proof->formula, out, std::forward<Printer>(printer)...) || !print('\n', out))
						return false;
				}
			}

			if (!print("  Elements: {", out)) return false;
			for (unsigned int k = 0; k < sets[i].element_count(); k++) {
				if (k != 0 && !print(", ", out)) return false;
				if (sets[i].arity > 1 && !print('(', out)) return false;
				for (unsigned int j = 0; j < sets[i].arity; j++) {
					if (j != 0 && !print(", ", out)) return false;
					if (!print(sets[i].elements[k * sets[i].arity + j], out, std::forward<Printer>(printer)...))
						return false;
				}
				if (sets[i].arity > 1 && !print(')', out)) return false;
			}
			if (!print("}\n", out)) return false;

			if (PrintProvableElements) {
				if (!print("  Provable elements: {", out)) return false;
				for (unsigned int k = 0; k < sets[i].provable_elements.length; k++) {
					if (k != 0 && !print(", ", out)) return false;
					if (!print(sets[i].provable_elements[k], out, std::forward<Printer>(printer)...))
						return false;
				}
				if (!print("}\n", out)) return false;
			}
		}
		return true;
	}

	bool check_freeable_sets() const {
		bool success = true;
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axioms.data == nullptr) continue;
			if (is_freeable(i)) {
				fprintf(stderr, "set_reasoning.check_freeable_sets: Set with ID %u is freeable.\n", i);
				success = false;
			}
		}
		return success;
	}

	bool are_set_sizes_valid() const {
		bool success = true;
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axioms.data == nullptr) continue;
			unsigned int min_set_size, max_set_size;
			if (!set_size_bounds(i, min_set_size, max_set_size))
				return false;
			if (sets[i].set_size < min_set_size || sets[i].set_size > max_set_size) {
				fprintf(stderr, "set_reasoning.are_set_sizes_valid WARNING: Set with ID %u has size (%u) outside the bounds computed by `set_reasoning.set_size_bounds` (min: %u, max: %u).\n", i, sets[i].set_size, min_set_size, max_set_size);
				success = false;
			}
		}
		return success;
	}

	bool check_set_ids() const {
		bool success = true;
		unsigned int actual_set_count = 0;
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axioms.data == nullptr) continue;

			bool contains;
			unsigned int set_id = set_ids.get(*sets[i].set_formula(), contains);
			if (!contains) {
				fprintf(stderr, "set_reasoning.check_set_ids WARNING: Set with ID %u has set formula '", i);
				print(*sets[i].set_formula(), stderr); fprintf(stderr, "' that does not exist in the `set_ids` map.\n");
				success = false;
			} else if (set_id != i) {
				fprintf(stderr, "set_reasoning.check_set_ids WARNING: Set with ID %u has set formula '", i);
				print(*sets[i].set_formula(), stderr); fprintf(stderr, "' that is mapped to %u in the `set_ids` map.\n", set_id);
				success = false;
			}
			actual_set_count++;
		}
		if (actual_set_count != set_ids.table.size) {
			fprintf(stderr, "set_reasoning.check_set_ids WARNING: `set_id` has size %u but there are a total of %u sets.\n", set_ids.table.size, actual_set_count);
			success = false;
		}

		for (const auto& entry : set_ids) {
			if (sets[entry.value].size_axioms.data == nullptr) {
				fprintf(stderr, "set_reasoning.check_set_ids WARNING: `set_ids` map contains an entry from formula '");
				print(entry.key, stderr);
				fprintf(stderr, "' to %u, but the set with ID %u doesn't exist (the size axiom is null).\n", entry.value, entry.value);
				success = false;
			} else if (*sets[entry.value].set_formula() != entry.key) {
				fprintf(stderr, "set_reasoning.check_set_ids WARNING: `set_ids` map contains an entry from formula '");
				print(entry.key, stderr); fprintf(stderr, "' to %u, but the set with ID %u has set formula '", entry.value, entry.value);
				print(*sets[entry.value].set_formula(), stderr); fprintf(stderr, "'.\n");
				success = false;
			}
		}
		return success;
	}

	bool check_symbols_in_formulas() const {
		hash_multiset<unsigned int> computed_symbols_in_formulas(256);
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axioms.data == nullptr) continue;

			array_multiset<unsigned int> symbols(16);
			Formula* formula = sets[i].set_formula();
			if (!get_constants(*formula, symbols)) return false;
			computed_symbols_in_formulas.add<true>(symbols);
		}

		bool success = true;
		for (const auto& entry : computed_symbols_in_formulas.counts) {
			bool contains;
			unsigned int count = symbols_in_formulas.counts.get(entry.key, contains);
			if (!contains) count = 0;
			if (entry.value != count) {
				fprintf(stderr, "set_reasoning.check_symbols_in_formulas WARNING: `symbols_in_formulas` for symbol `%u` has count %u, but only %u symbols appear in set formulas.\n", entry.key, count, entry.value);
				success = false;
			}
		} for (const auto& entry : symbols_in_formulas.counts) {
			if (computed_symbols_in_formulas.counts.table.contains(entry.key))
				continue;
			if (entry.value != 0) {
				fprintf(stderr, "set_reasoning.check_symbols_in_formulas WARNING: `symbols_in_formulas` for symbol `%u` has count %u, but only 0 symbols appear in set formulas.\n", entry.key, entry.value);
				success = false;
			}
		}
		return success;
	}

	bool are_provable_elements_valid() const {
		bool success = true;
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axioms.data == nullptr)
				continue;

			const tuple_element* elements_src = sets[i].elements.data;
			for (unsigned int k = 0; k < sets[i].element_count(); k++) {
				/* check for any duplicate elements */
				for (unsigned int l = k + 1; l < sets[i].element_count(); l++) {
					bool is_equal = true;
					for (uint_fast8_t m = 0; m < sets[i].arity; m++) {
						if (*(elements_src + (k * sets[i].arity) + m) != *(elements_src + (l * sets[i].arity) + m)) {
							is_equal = false;
							break;
						}
					}
					if (is_equal) {
						fprintf(stderr, "set_reasoning.are_provable_elements_valid WARNING: Set with ID %u has duplicate element ", i);
						if (sets[i].arity > 1 && !print('(', stderr)) return false;
						for (uint_fast8_t j = 0; j < sets[i].arity; j++) {
							if (j != 0 && !print(", ", stderr)) return false;
							if (!print(sets[i].elements[k * sets[i].arity + j], stderr, *debug_terminal_printer))
								return false;
						}
						if (sets[i].arity > 1 && !print(')', stderr)) return false;
						fprintf(stderr, ".\n");
						success = false;
						break;
					}
				}

				/* check if element k is in the list of provable_elements */
				bool is_element = false;
				for (const tuple& tup : sets[i].provable_elements) {
					bool is_equal = true;
					for (uint_fast8_t l = 0; l < sets[i].arity; l++) {
						if (*(elements_src + (k * sets[i].arity) + l) != tup.elements[l]) {
							is_equal = false;
							break;
						}
					}
					if (is_equal) {
						is_element = true;
						break;
					}
				}

				if (!is_element) {
					fprintf(stderr, "set_reasoning.are_provable_elements_valid WARNING: Set with ID %u has element ", i);
					if (sets[i].arity > 1 && !print('(', stderr)) return false;
					for (uint_fast8_t j = 0; j < sets[i].arity; j++) {
						if (j != 0 && !print(", ", stderr)) return false;
						if (!print(sets[i].elements[k * sets[i].arity + j], stderr, *debug_terminal_printer))
							return false;
					}
					if (sets[i].arity > 1 && !print(')', stderr)) return false;
					fprintf(stderr, " that is missing from `provable_elements`.\n");
					success = false;
				}
			}

			if (sets[i].element_count() > sets[i].provable_elements.length) {
				fprintf(stderr, "set_reasoning.are_provable_elements_valid WARNING: Set with ID %u has more elements than provable_elements.\n", i);
				success = false;
			} if (sets[i].element_count() > sets[i].provable_elements.capacity) {
				fprintf(stderr, "set_reasoning.are_provable_elements_valid WARNING: Set with ID %u has more elements than provable_elements has capacity.\n", i);
				success = false;
			}
		}
		return success;
	}
};

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool get_proof_map(
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		hash_map<const typename ProofCalculus::Proof*, unsigned int>& proof_map,
		hash_map<const typename ProofCalculus::Language*, unsigned int>& formula_map)
{
	for (unsigned int i = 1; i < sets.set_count + 1; i++) {
		if (sets.sets[i].size_axioms.data == nullptr) continue;
		if (!get_proof_map(sets.sets[i], proof_map, formula_map)
		 || !get_proof_map(sets.extensional_graph.vertices[i], proof_map, formula_map))
			return false;
	}

	for (const auto& entry : sets.set_ids)
		if (!get_formula_map(entry.key, formula_map)) return false;
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename Stream>
bool read(
		set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		Stream& in,
		typename ProofCalculus::Proof** proofs,
		typename ProofCalculus::Language** formulas)
{
	typedef typename ProofCalculus::Language Formula;

	unsigned int initialized_set_count;
	if (!read(initialized_set_count, in))
		return false;
	sets.capacity = ((size_t) 1) << (core::log2(initialized_set_count == 0 ? 1 : initialized_set_count) + 1);
	sets.intensional_graph.vertices = (intensional_set_vertex*) malloc(sizeof(intensional_set_vertex) * sets.capacity);
	if (sets.intensional_graph.vertices == nullptr) {
		fprintf(stderr, "read ERROR: Insufficient memory for `set_reasoning.intensional_set_graph.vertices`.\n");
		return false;
	}
	sets.extensional_graph.vertices = (extensional_set_vertex<ProofCalculus>*) malloc(sizeof(extensional_set_vertex<ProofCalculus>) * sets.capacity);
	if (sets.extensional_graph.vertices == nullptr) {
		fprintf(stderr, "read ERROR: Insufficient memory for `set_reasoning.extensional_graph.vertices`.\n");
		free(sets.intensional_graph.vertices); return false;
	}
	sets.sets = (set_info<BuiltInConstants, ProofCalculus>*) malloc(sizeof(set_info<BuiltInConstants, ProofCalculus>) * sets.capacity);
	if (sets.sets == nullptr) {
		fprintf(stderr, "read ERROR: Insufficient memory for `set_reasoning.extensional_graph.vertices`.\n");
		free(sets.intensional_graph.vertices); free(sets.extensional_graph.vertices);
		return false;
	}
	for (unsigned int i = 1; i < sets.capacity; i++)
		sets.sets[i].size_axioms.data = nullptr;
	sets.set_count = 0;

	for (unsigned int i = 0; i < initialized_set_count; i++) {
		unsigned int index;
		if (!read(index, in)
		 || !read(sets.sets[index], in, proofs))
		{
			for (unsigned int j = 1; j < sets.set_count + 1; j++) {
				free(sets.intensional_graph.vertices[j]);
				free(sets.extensional_graph.vertices[j]);
				free(sets.sets[j]);
			}
			free(sets.intensional_graph.vertices);
			free(sets.extensional_graph.vertices);
			free(sets.sets); return false;
		} else if (!read(sets.intensional_graph.vertices[index], in)) {
			free(sets.sets[index]);
			for (unsigned int j = 1; j < sets.set_count + 1; j++) {
				free(sets.intensional_graph.vertices[j]);
				free(sets.extensional_graph.vertices[j]);
				free(sets.sets[j]);
			}
			free(sets.intensional_graph.vertices);
			free(sets.extensional_graph.vertices);
			free(sets.sets); return false;
		} else if (!read(sets.extensional_graph.vertices[index], in, proofs)) {
			free(sets.intensional_graph.vertices[index]);
			free(sets.sets[index]);
			for (unsigned int j = 1; j < sets.set_count + 1; j++) {
				free(sets.intensional_graph.vertices[j]);
				free(sets.extensional_graph.vertices[j]);
				free(sets.sets[j]);
			}
			free(sets.intensional_graph.vertices);
			free(sets.extensional_graph.vertices);
			free(sets.sets); return false;
		}
		sets.set_count = index;
	}

	decltype(sets.set_ids.table.size) set_id_count;
	if (!read(set_id_count, in)
	 || !hash_map_init(sets.set_ids, 1 << (core::log2(RESIZE_THRESHOLD_INVERSE * (set_id_count == 0 ? 1 : set_id_count)) + 1)))
	{
		for (unsigned int j = 1; j < sets.set_count + 1; j++) {
			free(sets.intensional_graph.vertices[j]);
			free(sets.extensional_graph.vertices[j]);
			free(sets.sets[j]);
		}
		free(sets.intensional_graph.vertices);
		free(sets.extensional_graph.vertices);
		free(sets.sets); return false;
	}
	Formula& formula = *((Formula*) alloca(sizeof(Formula)));
	for (unsigned int i = 0; i < set_id_count; i++) {
		if (!read(formula, in, formulas)) {
			for (unsigned int j = 1; j < sets.set_count + 1; j++) {
				free(sets.intensional_graph.vertices[j]);
				free(sets.extensional_graph.vertices[j]);
				free(sets.sets[j]);
			} for (auto entry : sets.set_ids)
				free(entry.key);
			free(sets.intensional_graph.vertices);
			free(sets.extensional_graph.vertices);
			free(sets.sets); free(sets.set_ids);
			return false;
		}
		formula.reference_count = 1;
		unsigned int bucket = sets.set_ids.table.index_to_insert(formula);
		if (!read(sets.set_ids.values[bucket], in)) {
			free(formula);
			for (unsigned int j = 1; j < sets.set_count + 1; j++) {
				free(sets.intensional_graph.vertices[j]);
				free(sets.extensional_graph.vertices[j]);
				free(sets.sets[j]);
			} for (auto entry : sets.set_ids)
				free(entry.key);
			free(sets.intensional_graph.vertices);
			free(sets.extensional_graph.vertices);
			free(sets.sets); free(sets.set_ids);
			return false;
		}
		move(formula, sets.set_ids.table.keys[bucket]);
		sets.set_ids.table.size++;
	}

	if (!read(sets.symbols_in_formulas, in)) {
		for (unsigned int j = 1; j < sets.set_count + 1; j++) {
			free(sets.intensional_graph.vertices[j]);
			free(sets.extensional_graph.vertices[j]);
			free(sets.sets[j]);
		} for (auto entry : sets.set_ids)
			free(entry.key);
		free(sets.intensional_graph.vertices);
		free(sets.extensional_graph.vertices);
		free(sets.sets); free(sets.set_ids);
		return false;
	}
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename Stream>
bool write(
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets, Stream& out,
		const hash_map<const typename ProofCalculus::Proof*, unsigned int>& proof_map,
		const hash_map<const typename ProofCalculus::Language*, unsigned int>& formula_map)
{
	unsigned int initialized_set_count = 0;
	for (unsigned int i = 1; i < sets.set_count + 1; i++)
		if (sets.sets[i].size_axioms.data != nullptr) initialized_set_count++;
	if (!write(initialized_set_count, out))
		return false;
	for (unsigned int i = 1; i < sets.set_count + 1; i++) {
		if (sets.sets[i].size_axioms.data == nullptr) continue;
		if (!write(i, out)
		 || !write(sets.sets[i], out, proof_map)
		 || !write(sets.intensional_graph.vertices[i], out)
		 || !write(sets.extensional_graph.vertices[i], out, proof_map))
			return false;
	}

	if (!write(sets.set_ids.table.size, out))
		return false;
	for (const auto& entry : sets.set_ids) {
		if (!write(entry.key, out, formula_map)
		 || !write(entry.value, out))
			return false;
	}

	return write(sets.symbols_in_formulas, out);
}

struct clique_search_state {
	unsigned int* clique;
	unsigned int clique_count;
	unsigned int* neighborhood;
	unsigned int neighborhood_count;
	unsigned int* X;
	unsigned int X_count;
	unsigned int next_set;
	int priority;

	inline int get_priority() const {
		return priority;
	}

	static inline void free(clique_search_state& state) {
		core::free(state.clique);
		core::free(state.neighborhood);
		core::free(state.X);
	}
};

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool init(clique_search_state& state,
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		const unsigned int* clique, unsigned int clique_count,
		const unsigned int* neighborhood, unsigned int neighborhood_count,
		const unsigned int* X, unsigned int X_count,
		unsigned int next_set_to_expand, int priority)
{
	state.clique = (unsigned int*) malloc(max((size_t) 1, sizeof(unsigned int) * clique_count));
	if (state.clique == NULL) {
		fprintf(stderr, "init ERROR: Insufficient memory for `clique_search_state.clique`.");
		return false;
	}
	for (unsigned int i = 0; i < clique_count; i++)
		state.clique[i] = clique[i];
	state.clique_count = clique_count;

	state.neighborhood = (unsigned int*) malloc(max((size_t) 1, sizeof(unsigned int) * neighborhood_count));
	if (state.neighborhood == NULL) {
		fprintf(stderr, "init ERROR: Insufficient memory for `clique_search_state.neighborhood`.");
		core::free(state.clique); return false;
	}
	for (unsigned int i = 0; i < neighborhood_count; i++)
		state.neighborhood[i] = neighborhood[i];
	state.neighborhood_count = neighborhood_count;

	state.X = (unsigned int*) malloc(max((size_t) 1, sizeof(unsigned int) * X_count));
	if (state.X == NULL) {
		fprintf(stderr, "init ERROR: Insufficient memory for `clique_search_state.X`.");
		core::free(state.clique); core::free(state.neighborhood); return false;
	}
	for (unsigned int i = 0; i < X_count; i++)
		state.X[i] = X[i];
	state.X_count = X_count;

	state.next_set = next_set_to_expand;
	state.priority = priority;
	return true;
}

template<typename Stream>
bool print(const clique_search_state& state, Stream& out) {
	return print("  clique: ", out) && print<unsigned int, '{', '}'>(state.clique, state.clique_count, out) && print('\n', out)
		&& print("  neighborhood: ", out) && print<unsigned int, '{', '}'>(state.neighborhood, state.neighborhood_count, out) && print('\n', out)
		&& print("  X: ", out) && print<unsigned int, '{', '}'>(state.X, state.X_count, out) && print('\n', out)
		&& print("  next_set: ", out) && print(state.next_set, out) && print('\n', out);
}

struct ancestor_clique_search_state {
	clique_search_state clique_state;
	unsigned int ancestor;

	inline int get_priority() const {
		return clique_state.priority;
	}

	static inline void free(ancestor_clique_search_state& state) {
		core::free(state.clique_state);
	}
};

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool init(ancestor_clique_search_state& state,
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		const unsigned int* clique, unsigned int clique_count,
		const unsigned int* neighborhood, unsigned int neighborhood_count,
		const unsigned int* X, unsigned int X_count,
		unsigned int next_set_to_expand, int neighborhood_size,
		unsigned int ancestor)
{
	if (!init(state.clique_state, sets, clique, clique_count, neighborhood, neighborhood_count,
			X, X_count, next_set_to_expand, neighborhood_size - sets.sets[ancestor].set_size))
		return false;
	state.ancestor = ancestor;
	return true;
}

template<typename StateData>
struct search_state {
	StateData* state;

	static inline void free(search_state<StateData>& state) {
		core::free(*state.state);
		core::free(state.state);
	}
};

template<typename StateData>
inline bool operator < (const search_state<StateData>& first, const search_state<StateData>& second) {
	return first.state->get_priority() < second.state->get_priority();
}

template<typename StateData, typename Stream>
inline bool print(const search_state<StateData>& state, Stream& out) {
	return print(state.state, out)
		&& print("  priority: ", out) && print(state.priority, out) && print('\n', out);
}

template<typename StateData>
struct search_queue {
	std::multiset<search_state<StateData>> queue;
	int last_priority;

	search_queue() : last_priority(INT_MAX) { }

	~search_queue() {
		for (search_state<StateData> state : queue)
			core::free(state);
	}

	inline bool is_empty() const {
		return queue.empty();
	}

	inline unsigned int size() const {
		return queue.size();
	}

	inline void push(const search_state<StateData>& state) {
#if !defined(NDEBUG)
		if (state.state->get_priority() > last_priority)
			fprintf(stderr, "search_queue.push WARNING: Search is not monotonic.\n");
#endif
		queue.insert(state);
	}

	inline search_state<StateData> pop(unsigned int iteration) {
		auto last = queue.cend(); last--;
		search_state<StateData> state = *last;
		queue.erase(last);

		last_priority = state.state->get_priority();
		return state;
	}

	inline int priority() const {
		auto last = queue.cend(); last--;
		search_state<StateData> state = *last;
		return state.state->get_priority();
	}
};

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool has_descendant(
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int root, unsigned int vertex)
{
	return sets.sets[root].descendants.contains(vertex);
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
bool expand_clique_search_state(
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		array<unsigned int>& neighborhood,
		const unsigned int* X, const unsigned int X_count,
		hash_set<unsigned int>& visited,
		unsigned int set_to_expand,
		unsigned int root)
{
	if (!visited.add(root)) return false;

	if (sets.are_disjoint(set_to_expand, root)) {
		for (unsigned int i = 0; i < neighborhood.length; i++)
			if (has_descendant(sets, neighborhood[i], root)) { neighborhood.remove(i); i--; }

		if (!neighborhood.ensure_capacity(neighborhood.length + 1))
			return false;

#if !defined(NDEBUG)
		if (neighborhood.contains(root))
			fprintf(stderr, "expand_clique_search_state WARNING: `neighborhood` contains `root`.\n");
#endif

		neighborhood[neighborhood.length++] = root;
	} else {
		for (const auto& entry : sets.extensional_graph.vertices[root].children) {
			if (sets.sets[entry.key].set_size == 0) continue;
			bool was_visited = visited.contains(entry.key);
			for (unsigned int i = 0; !was_visited && i < X_count; i++)
				if (has_descendant(sets, X[i], entry.key)) was_visited = true;
			for (unsigned int neighbor : neighborhood)
				if (has_descendant(sets, neighbor, entry.key)) was_visited = true;
			if (was_visited) continue;
			if (!expand_clique_search_state(sets, neighborhood, X, X_count, visited, set_to_expand, entry.key)) return false;
		} for (unsigned int child : sets.intensional_graph.vertices[root].children) {
			if (sets.sets[child].set_size == 0) continue;
			bool was_visited = visited.contains(child);
			for (unsigned int i = 0; !was_visited && i < X_count; i++)
				if (has_descendant(sets, X[i], child)) was_visited = true;
			for (unsigned int neighbor : neighborhood)
				if (has_descendant(sets, neighbor, child)) was_visited = true;
			if (was_visited) continue;
			if (!expand_clique_search_state(sets, neighborhood, X, X_count, visited, set_to_expand, child)) return false;
		}
	}
	return true;
}

template<bool RecurseChildren, bool TestCompletion, bool ReturnOnCompletion, typename StateData,
	typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename... StateArgs>
bool process_clique_search_state(
		search_queue<StateData>& queue,
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		const clique_search_state& state,
		unsigned int*& clique, unsigned int& clique_count,
		int min_priority, StateArgs&&... state_args)
{
	unsigned int* new_clique = (unsigned int*) malloc(sizeof(unsigned int) * (state.clique_count + 1));
	if (new_clique == NULL) {
		fprintf(stderr, "process_clique_search_state ERROR: Insufficient memory for `new_clique`.\n");
		return false;
	}
	for (unsigned int i = 0; i < state.clique_count; i++)
		new_clique[i] = state.clique[i];
	new_clique[state.clique_count] = state.next_set;

	array<unsigned int> neighborhood(16);
	hash_set<unsigned int> visited(16);
	for (unsigned int i = 0; i < state.X_count; i++) {
		if (!expand_clique_search_state(sets, neighborhood, neighborhood.data, neighborhood.length, visited, state.next_set, state.X[i])) {
			free(new_clique);
			return false;
		}
	}
	unsigned int new_X_count = neighborhood.length;
	for (unsigned int i = 0; i < state.neighborhood_count; i++) {
		if (!expand_clique_search_state(sets, neighborhood, neighborhood.data, new_X_count, visited, state.next_set, state.neighborhood[i])) {
			free(new_clique);
			return false;
		}
	}

	int priority = sets.sets[state.next_set].set_size;
	for (unsigned int i = 0; i < state.clique_count; i++)
		priority += sets.sets[state.clique[i]].set_size;
	for (unsigned int i = new_X_count; i < neighborhood.length; i++)
		priority += sets.sets[neighborhood[i]].set_size;

	for (unsigned int i = new_X_count; i < neighborhood.length; i++) {
		if (priority < min_priority) break;
		search_state<StateData> new_state;
		new_state.state = (StateData*) malloc(sizeof(StateData));
		if (new_state.state == NULL) {
			free(new_clique);
			return false;
		} else if (!init(*new_state.state, sets, new_clique, state.clique_count + 1,
				neighborhood.data + i + 1, neighborhood.length - i - 1, neighborhood.data, i, neighborhood[i],
				min(state.priority, priority), std::forward<StateArgs>(state_args)...))
		{
			free(new_state.state); free(new_clique);
			return false;
		}

		queue.push(new_state);
		priority -= sets.sets[neighborhood[i]].set_size;
	}

	if (TestCompletion && neighborhood.length == 0 && state.X_count == 0) {
		/* `expand_clique_search_state` did not add any states to the queue, meaning this clique is maximal */
		clique = new_clique;
		clique_count = state.clique_count + 1;
		if (ReturnOnCompletion) return true;
	} else {
		free(new_clique);
	}

	if (RecurseChildren) {
		unsigned int old_neighborhood_length = neighborhood.length;
		const auto& extensional_children = sets.extensional_graph.vertices[state.next_set].children;
		const array<unsigned int>& intensional_children = sets.intensional_graph.vertices[state.next_set].children;
		array<unsigned int> children(max((size_t) 1, extensional_children.size + intensional_children.length));
		for (const auto& entry : extensional_children) {
			/* make sure we don't recurse into an ancestor of the current node */
			if (sets.sets[entry.key].set_size == 0) continue;
			children[children.length++] = entry.key;
		} for (unsigned int i = 0; i < intensional_children.length; i++) {
			if (extensional_children.contains(intensional_children[i]) || sets.sets[intensional_children[i]].set_size == 0) continue;
			children[children.length++] = intensional_children[i];
		}

		if (!neighborhood.ensure_capacity(neighborhood.length + children.length))
			return false;

		priority -= sets.sets[state.next_set].set_size;
		int priority_without_parent = priority;

		for (unsigned int i = 0; i < children.length; i++) {
			unsigned int child = children[i];
			neighborhood.length = old_neighborhood_length;
			priority = priority_without_parent + sets.sets[child].set_size;

			/* only add the children to neighborhood that are disjoint with `child` */
			for (unsigned int j = 0; j < i; j++) {
				if (!sets.are_disjoint(child, children[j])) continue;
				neighborhood[neighborhood.length++] = children[j];
				priority += sets.sets[children[j]].set_size;
			}

			search_state<StateData> new_state;
			new_state.state = (StateData*) malloc(sizeof(StateData));
			if (new_state.state == NULL) {
				return false;
			} else if (!init(*new_state.state, sets, state.clique, state.clique_count,
					neighborhood.data + old_neighborhood_length,
					neighborhood.length - old_neighborhood_length,
					neighborhood.data, old_neighborhood_length, child,
					min(state.priority, priority), std::forward<StateArgs>(state_args)...))
			{
				free(new_state.state);
				return false;
			}
			queue.push(new_state);
		}
	}

	return true;
}

/* this is the Bron-Kerbosch algorithm */
template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
bool find_largest_disjoint_subset_clique(
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int root, unsigned int*& clique, unsigned int& clique_count,
		int min_priority)
{
	search_queue<clique_search_state> queue;

	clique_search_state initial_state;
	if (!init(initial_state, sets, NULL, 0, NULL, 0, NULL, 0, root, INT_MAX)) {
		return false;
	} else if (!process_clique_search_state<true, false, false>(queue, sets, initial_state, clique, clique_count, min_priority)) {
		free(initial_state);
		return false;
	}
	free(initial_state);

	bool success = true;
	clique = NULL; clique_count = 0;
	int best_clique_score = INT_MIN;
	for (unsigned int iteration = 0; success && !queue.is_empty() && queue.priority() > best_clique_score; iteration++) {
		search_state<clique_search_state> next = queue.pop(iteration);
		unsigned int* completed_clique = NULL; unsigned int completed_clique_count;
		success = process_clique_search_state<true, true, true>(queue, sets, *next.state, completed_clique, completed_clique_count, min_priority);
		free(next);

		if (completed_clique != NULL) {
			int priority = 0;
			for (unsigned int i = 0; i < completed_clique_count; i++)
				priority += sets.sets[completed_clique[i]].set_size;
			if (priority > best_clique_score && priority >= min_priority) {
				if (clique != NULL) free(clique);
				clique = completed_clique;
				clique_count = completed_clique_count;
				best_clique_score = priority;
			} else {
				free(completed_clique);
			}
		}
	}

#if !defined(NDEBUG)
	for (unsigned int i = 0; clique != NULL && i < clique_count; i++)
		if (sets.sets[clique[i]].set_size == 0)
			fprintf(stderr, "find_largest_disjoint_subset_clique WARNING: `clique` contains an empty set.\n");
#endif

	return success;
}

inline bool add_non_ancestor_neighbor(
		unsigned int child, bool& changed, array<unsigned int>& neighborhood,
		hash_map<unsigned int, array<unsigned int>>& non_ancestor_neighborhood)
{
	bool contains;
	array<unsigned int>& child_neighborhood = non_ancestor_neighborhood.get(child, contains);
	if (contains) {
		for (unsigned int child_neighbor : child_neighborhood) {
			if (neighborhood.contains(child_neighbor)) continue;
			if (!neighborhood.add(child_neighbor)) return false;
			changed = true;
		}
	} else {
		if (neighborhood.contains(child)) return true;
		if (!neighborhood.add(child)) return false;
		changed = true;
	}
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool add_non_ancestor_neighbors(
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int node, bool& changed,
		hash_map<unsigned int, array<unsigned int>>& non_ancestor_neighborhood)
{
	array<unsigned int>& neighborhood = non_ancestor_neighborhood.get(node);
	for (const auto& entry : sets.extensional_graph.vertices[node].children)
		if (sets.sets[entry.key].set_size > 0 && !add_non_ancestor_neighbor(entry.key, changed, neighborhood, non_ancestor_neighborhood)) return false;
	for (unsigned int child : sets.intensional_graph.vertices[node].children)
		if (sets.sets[child].set_size > 0 && !add_non_ancestor_neighbor(child, changed, neighborhood, non_ancestor_neighborhood)) return false;
	return true;
}

struct default_parents {
	template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename Function>
	bool for_each_parent(const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets, unsigned int set, Function apply) const {
		for (const auto& entry : sets.extensional_graph.vertices[set].parents)
			if (!apply(entry.key)) return false;
		for (unsigned int parent : sets.intensional_graph.vertices[set].parents)
			if (!apply(parent)) return false;
		return true;
	}

	template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
	inline unsigned int count(const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets, unsigned int set) const {
		return sets.extensional_graph.vertices[set].parents.size
			 + sets.intensional_graph.vertices[set].parents.length;
	}
};

struct singleton_parent {
	unsigned int parent;

	template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename Function>
	bool for_each_parent(const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets, unsigned int set, Function apply) const {
		return apply(parent);
	}

	template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
	constexpr unsigned int count(const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets, unsigned int set) const {
		return 1;
	}
};

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename SetParents>
bool get_non_ancestor_neighborhoods(
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int set, const SetParents& set_parents,
		hash_map<unsigned int, array<unsigned int>>& non_ancestor_neighborhood)
{
	/* first collect all ancestors of `set` */
	unsigned int parent_count = set_parents.count(sets, set);
	array<unsigned int> stack((parent_count == 0) ? 1 : (1 << (core::log2(parent_count) + 1)));
	auto init_non_ancestor_neighborhood = [&](unsigned int current) {
		if (!non_ancestor_neighborhood.check_size())
			return false;

		bool contains; unsigned int bucket;
		array<unsigned int>& siblings = non_ancestor_neighborhood.get(current, contains, bucket);
		if (!contains) {
			if (!array_init(siblings, 4)) return false;
			non_ancestor_neighborhood.table.keys[bucket] = current;
			non_ancestor_neighborhood.table.size++;
		}
		return true;
	};
	if (!init_non_ancestor_neighborhood(set)) {
		for (auto entry : non_ancestor_neighborhood) free(entry.value);
		return false;
	}
	stack[stack.length++] = set;
	while (stack.length > 0) {
		unsigned int current = stack.pop();

		for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
			if (entry.key != set && non_ancestor_neighborhood.table.contains(entry.key)) continue;
			if (!init_non_ancestor_neighborhood(entry.key) || !stack.add(entry.key)) {
				for (auto entry : non_ancestor_neighborhood) free(entry.value);
				return false;
			}
		} for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
			if (parent != set && non_ancestor_neighborhood.table.contains(parent)) continue;
			if (!init_non_ancestor_neighborhood(parent) || !stack.add(parent)) {
				for (auto entry : non_ancestor_neighborhood) free(entry.value);
				return false;
			}
		}
	}

	auto process_ancestor_parent = [&](unsigned int parent) {
		bool parent_changed = false;
		if (!add_non_ancestor_neighbors(sets, parent, parent_changed, non_ancestor_neighborhood))
			return false;
		if (parent_changed && !stack.add(parent)) return false;
		return true;
	};
	if (!set_parents.for_each_parent(sets, set, process_ancestor_parent)) {
		for (auto entry : non_ancestor_neighborhood) free(entry.value);
		return false;
	}
	while (stack.length > 0) {
		unsigned int current = stack.pop();
		for (const auto& entry : sets.extensional_graph.vertices[current].parents) {
			if (!process_ancestor_parent(entry.key)) {
				for (auto entry : non_ancestor_neighborhood) free(entry.value);
				return false;
			}
		} for (unsigned int parent : sets.intensional_graph.vertices[current].parents) {
			if (!process_ancestor_parent(parent)) {
				for (auto entry : non_ancestor_neighborhood) free(entry.value);
				return false;
			}
		}
	}

	bool contains; unsigned int bucket;
	free(non_ancestor_neighborhood.get(set, contains, bucket));
	non_ancestor_neighborhood.remove_at(bucket);
	return true;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename SetParents>
bool find_largest_disjoint_clique_with_set(
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int set, const SetParents& set_parents,
		unsigned int*& clique, unsigned int& clique_count,
		unsigned int& ancestor_of_clique, int min_priority)
{
	hash_map<unsigned int, array<unsigned int>> non_ancestor_neighborhood(32);
	if (!get_non_ancestor_neighborhoods(sets, set, set_parents, non_ancestor_neighborhood))
		return false;

	search_queue<ancestor_clique_search_state> queue;
	clique = NULL; clique_count = 0;
	int best_clique_score = INT_MIN;
	for (const auto& entry : non_ancestor_neighborhood) {
		ancestor_clique_search_state initial_state;
		initial_state.ancestor = entry.key;
		unsigned int* completed_clique = NULL; unsigned int completed_clique_count;
		if (!init(initial_state.clique_state, sets, NULL, 0, entry.value.data, entry.value.length, NULL, 0, set, INT_MAX)) {
			for (auto entry : non_ancestor_neighborhood) free(entry.value);
			return false;
		} else if (!process_clique_search_state<false, true, false>(
				queue, sets, initial_state.clique_state, completed_clique,
				completed_clique_count, min_priority, initial_state.ancestor))
		{
			free(initial_state);
			for (auto entry : non_ancestor_neighborhood) free(entry.value);
			return false;
		}
		free(initial_state);

		if (completed_clique != NULL) {
			/* this only happens if `set` is disjoint with everything in the initial neighborhood */
			int priority = sets.sets[set].set_size - sets.sets[entry.key].set_size;
			if (priority > best_clique_score && priority >= min_priority) {
				if (clique != NULL) free(clique);
				clique = completed_clique;
				clique_count = completed_clique_count;
				ancestor_of_clique = entry.key;
				best_clique_score = priority;
			} else {
				free(completed_clique);
			}
		}
	}
	for (auto entry : non_ancestor_neighborhood) free(entry.value);

	bool success = true;
	for (unsigned int iteration = 0; success && !queue.is_empty() && queue.priority() > best_clique_score; iteration++) {
		search_state<ancestor_clique_search_state> next = queue.pop(iteration);
		unsigned int* completed_clique = NULL; unsigned int completed_clique_count;
		success = process_clique_search_state<true, true, true>(
				queue, sets, next.state->clique_state, completed_clique,
				completed_clique_count, min_priority, next.state->ancestor);

		if (completed_clique != NULL) {
			int priority = -((int) sets.sets[next.state->ancestor].set_size);
			for (unsigned int i = 0; i < completed_clique_count; i++)
				priority += sets.sets[completed_clique[i]].set_size;
			if (priority > best_clique_score && priority >= min_priority) {
				if (clique != NULL) free(clique);
				clique = completed_clique;
				clique_count = completed_clique_count;
				ancestor_of_clique = next.state->ancestor;
				best_clique_score = priority;
			} else {
				free(completed_clique);
			}
		}
		free(next);
	}

#if !defined(NDEBUG)
	for (unsigned int i = 0; clique != NULL && i < clique_count; i++)
		if (sets.sets[clique[i]].set_size == 0 && clique[i] != set)
			fprintf(stderr, "find_largest_disjoint_clique_with_set WARNING: `clique` contains an empty set.\n");
#endif

	return success;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool find_largest_disjoint_clique_with_set(
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int set, unsigned int*& clique, unsigned int& clique_count,
		unsigned int& ancestor_of_clique, int min_priority)
{
	default_parents parents;
	return find_largest_disjoint_clique_with_set(sets, set, parents, clique, clique_count, ancestor_of_clique, min_priority);
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool find_largest_disjoint_clique_with_set(
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int set, unsigned int parent,
		unsigned int*& clique, unsigned int& clique_count,
		unsigned int& ancestor_of_clique, int min_priority)
{
	singleton_parent set_parent = { parent };
	return find_largest_disjoint_clique_with_set(sets, set, set_parent, clique, clique_count, ancestor_of_clique, min_priority);
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer, typename FirstSetParents, typename SecondSetParents>
bool find_largest_disjoint_clique_with_edge(
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int first, unsigned int second,
		const FirstSetParents& first_set_parents,
		const SecondSetParents& second_set_parents,
		unsigned int*& clique, unsigned int& clique_count,
		unsigned int& ancestor_of_clique, int min_priority,
		int stop_priority)
{
	hash_map<unsigned int, array<unsigned int>> first_non_ancestor_neighborhood(32);
	hash_map<unsigned int, array<unsigned int>> second_non_ancestor_neighborhood(32);
	if (!get_non_ancestor_neighborhoods(sets, first, first_set_parents, first_non_ancestor_neighborhood)) {
		return false;
	} else if (!get_non_ancestor_neighborhoods(sets, second, second_set_parents, second_non_ancestor_neighborhood)) {
		for (auto entry : first_non_ancestor_neighborhood) free(entry.value);
		return false;
	}

	/* intersect the two non-ancestor neighborhoods */
	hash_map<unsigned int, array<unsigned int>> non_ancestor_neighborhood(RESIZE_THRESHOLD_INVERSE * (max(first_non_ancestor_neighborhood.table.size, second_non_ancestor_neighborhood.table.size) + 1));
	for (auto entry : first_non_ancestor_neighborhood) {
		bool contains;
		array<unsigned int>& first_neighborhood = entry.value;
		if (!first_neighborhood.contains(second))
			continue;
		array<unsigned int>& second_neighborhood = second_non_ancestor_neighborhood.get(entry.key, contains);
		if (contains) {
			unsigned int bucket;
			array<unsigned int>& intersection = non_ancestor_neighborhood.get(entry.key, contains, bucket);
			if (!array_init(intersection, max(first_neighborhood.length, second_neighborhood.length) + 1)) {
				for (auto entry : first_non_ancestor_neighborhood) free(entry.value);
				for (auto entry : second_non_ancestor_neighborhood) free(entry.value);
				for (auto entry : non_ancestor_neighborhood) free(entry.value);
				return false;
			}
			non_ancestor_neighborhood.table.keys[bucket] = entry.key;
			non_ancestor_neighborhood.table.size++;

			if (first_neighborhood.length > 1)
				insertion_sort(first_neighborhood);
			if (second_neighborhood.length > 1)
				insertion_sort(second_neighborhood);
			set_intersect(intersection, first_neighborhood, second_neighborhood);
			if (intersection.length != 0)
				move(intersection[0], intersection[intersection.length]);
			intersection[0] = second;
			intersection.length++;
		}
	}
	for (auto entry : first_non_ancestor_neighborhood) free(entry.value);
	for (auto entry : second_non_ancestor_neighborhood) free(entry.value);

	search_queue<ancestor_clique_search_state> queue;
	clique = NULL; clique_count = 0;
	for (const auto& entry : non_ancestor_neighborhood) {
		/* create search states that contain both `first` and `second` in the clique */
		unsigned int* new_clique = (unsigned int*) malloc(sizeof(unsigned int));
		if (new_clique == NULL) {
			fprintf(stderr, "find_largest_disjoint_clique_with_edge ERROR: Insufficient memory for `new_clique`.\n");
			for (auto entry : non_ancestor_neighborhood) free(entry.value);
			return false;
		}
		new_clique[0] = first;

		array<unsigned int> neighborhood(16);
		hash_set<unsigned int> visited(16);
		unsigned int new_X_count = neighborhood.length;
		for (unsigned int i = 0; i < entry.value.length; i++) {
			if (!expand_clique_search_state(sets, neighborhood, neighborhood.data, new_X_count, visited, first, entry.value[i])) {
				free(new_clique);
				return false;
			}
		}

		int priority = sets.sets[first].set_size;
		for (unsigned int i = new_X_count; i < neighborhood.length; i++)
			priority += sets.sets[neighborhood[i]].set_size;

		if (priority >= min_priority) {
			search_state<ancestor_clique_search_state> new_state;
			new_state.state = (ancestor_clique_search_state*) malloc(sizeof(ancestor_clique_search_state));
			priority += sets.sets[second].set_size;
			if (new_state.state == NULL) {
				for (auto entry : non_ancestor_neighborhood) free(entry.value);
				free(new_clique); return false;
			} else if (!init(*new_state.state, sets, new_clique, 1,
					neighborhood.data + 1, neighborhood.length - 1, neighborhood.data, 0,
					second, priority, entry.key))
			{
				for (auto entry : non_ancestor_neighborhood) free(entry.value);
				free(new_state.state); free(new_clique); return false;
			}
			queue.push(new_state);
		}
		free(new_clique);
	}
	for (auto entry : non_ancestor_neighborhood) free(entry.value);

	bool success = true;
	for (unsigned int iteration = 0; success && !queue.is_empty(); iteration++) {
		search_state<ancestor_clique_search_state> next = queue.pop(iteration);

		int priority = sets.sets[next.state->ancestor].set_size - sets.sets[next.state->clique_state.next_set].set_size;
		for (unsigned int i = 0; i < next.state->clique_state.clique_count; i++)
			priority -= sets.sets[next.state->clique_state.clique[i]].set_size;
		if (priority < stop_priority) {
			if (clique != NULL) free(clique);
			clique = (unsigned int*) malloc(sizeof(unsigned int) * (next.state->clique_state.clique_count + 1));
			memcpy(clique, next.state->clique_state.clique, sizeof(unsigned int) * next.state->clique_state.clique_count);
			clique[next.state->clique_state.clique_count] = next.state->clique_state.next_set;
			clique_count = next.state->clique_state.clique_count + 1;
			ancestor_of_clique = next.state->ancestor;
			free(next); break;
		}

		unsigned int* completed_clique = NULL; unsigned int completed_clique_count;
		success = process_clique_search_state<true, true, true>(
				queue, sets, next.state->clique_state, completed_clique,
				completed_clique_count, min_priority, next.state->ancestor);

		if (completed_clique != NULL)
			free(completed_clique);
		free(next);
	}

#if !defined(NDEBUG)
	for (unsigned int i = 0; clique != NULL && i < clique_count; i++)
		if (sets.sets[clique[i]].set_size == 0 && clique[i] != first && clique[i] != second)
			fprintf(stderr, "find_largest_disjoint_clique_with_set WARNING: `clique` contains an empty set.\n");
#endif

	return success;
}

template<typename BuiltInConstants, typename ProofCalculus, typename Canonicalizer>
inline bool find_largest_disjoint_clique_with_edge(
		const set_reasoning<BuiltInConstants, ProofCalculus, Canonicalizer>& sets,
		unsigned int first, unsigned int second,
		unsigned int*& clique, unsigned int& clique_count,
		unsigned int& ancestor_of_clique, int min_priority,
		int stop_priority)
{
	default_parents parents;
	return find_largest_disjoint_clique_with_edge(sets, first, second, parents, parents, clique, clique_count, ancestor_of_clique, min_priority, stop_priority);
}

#endif /* SET_GRAPH_H_ */
