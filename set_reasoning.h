#ifndef SET_GRAPH_H_
#define SET_GRAPH_H_

#include <core/array.h>
#include <core/map.h>
#include <stdint.h>
#include <set>

using namespace core;


/* forward declarations */

template<typename BuiltInConstants, typename Formula, typename ProofCalculus> struct set_reasoning;

template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
bool find_largest_disjoint_clique_with_set(
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>&,
		unsigned int, unsigned int*&, unsigned int&, unsigned int&, int);
template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
bool find_largest_disjoint_clique_with_set(
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>&,
		unsigned int, unsigned int, unsigned int*&, unsigned int&, unsigned int&, int);

template<typename T, typename Stream>
bool print(const hash_set<T>& set, Stream& out) {
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

template<typename ProofCalculus>
struct extensional_set_vertex
{
	typedef typename ProofCalculus::Proof Proof;

	array_map<unsigned int, Proof*> parents;
	array_map<unsigned int, Proof*> children;

	static inline void free(extensional_set_vertex<ProofCalculus>& vertex) {
		for (auto entry : vertex.children) {
			core::free(*entry.value); if (entry.value->reference_count == 0) core::free(entry.value);
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

struct intensional_set_graph {
	intensional_set_vertex* vertices;

	intensional_set_graph(unsigned int initial_capacity) {
		vertices = (intensional_set_vertex*) malloc(sizeof(intensional_set_vertex) * initial_capacity);
		if (vertices == NULL) exit(EXIT_FAILURE);
	}

	~intensional_set_graph() {
		free(vertices);
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
		free(vertices);
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
			for (const auto& entry : vertices[vertex_id].parents)
				remove_edge(entry.key, vertex_id);
			for (const auto& entry : vertices[vertex_id].children)
				remove_edge(vertex_id, entry.key);
		}
		core::free(vertices[vertex_id]);
	}

	inline bool add_edge(unsigned int parent, unsigned int child, Proof* axiom)
	{
		if (!vertices[parent].children.ensure_capacity(vertices[parent].children.size + 1)
		 || !vertices[child].parents.ensure_capacity(vertices[child].parents.size + 1))
			return false;
#if !defined(NDEBUG)
		if (vertices[parent].children.contains(child)) {
			fprintf(stderr, "extensional_edge.add_edge WARNING: The given edge already exists.\n");
			return false;
		}
#endif

		unsigned int index = vertices[parent].children.size;
		vertices[parent].children.keys[index] = child;
		vertices[parent].children.values[index] = axiom;
		vertices[parent].children.size++;

		index = vertices[child].parents.size;
		vertices[child].parents.keys[index] = parent;
		vertices[child].parents.values[index] = axiom;
		vertices[child].parents.size++;
		axiom->reference_count++;
		return true;
	}

	template<typename Formula>
	inline Proof* get_edge(unsigned int parent, unsigned int child,
			Formula* parent_formula, Formula* child_formula, bool& new_edge)
	{
		if (!vertices[parent].children.ensure_capacity(vertices[parent].children.size + 1)
		 || !vertices[child].parents.ensure_capacity(vertices[child].parents.size + 1))
			return NULL;
		unsigned int index = vertices[parent].children.index_of(child);
		if (index == vertices[parent].children.size) {
			new_edge = true;
			Formula* formula = Formula::new_for_all(1, Formula::new_if_then(child_formula, parent_formula));
			if (formula == NULL) return NULL;
			Proof* new_axiom = ProofCalculus::new_axiom(formula);
			free(*formula); if (formula->reference_count == 0) free(formula);
			if (new_axiom == NULL) return NULL;
			child_formula->reference_count++;
			parent_formula->reference_count++;

			vertices[parent].children.keys[index] = child;
			vertices[parent].children.values[index] = new_axiom;
			vertices[parent].children.size++;

			index = vertices[child].parents.size;
			vertices[child].parents.keys[index] = parent;
			vertices[child].parents.values[index] = new_axiom;
			vertices[child].parents.size++;
			new_axiom->reference_count++;
			return new_axiom;
		} else {
			new_edge = false;
			return vertices[parent].children.values[index];
		}
	}

	inline void remove_edge(unsigned int parent, unsigned int child) {
		unsigned int index = vertices[parent].children.index_of(child);
#if !defined(NDEBUG)
		if (index == vertices[parent].children.size)
			fprintf(stderr, "extensional_set_graph.remove_edge WARNING: The set %u is not in `vertices[%u].children`.\n", child, parent);
#endif
		Proof* axiom = vertices[parent].children.values[index];
		core::free(*axiom);
		if (axiom->reference_count == 0) core::free(axiom);
		vertices[parent].children.remove_at(index);
		index = vertices[child].parents.index_of(parent);
#if !defined(NDEBUG)
		if (index == vertices[child].parents.size)
			fprintf(stderr, "extensional_set_graph.remove_edge WARNING: The set %u is not in `vertices[%u].children`.\n", parent, child);
#endif
		vertices[child].parents.remove_at(index);
	}
};

bool get_ancestors(
		const intensional_set_graph& graph, unsigned int vertex,
		array_map<unsigned int, unsigned int>& ancestors)
{
	if (!ancestors.ensure_capacity(ancestors.size + 1)) return false;
	unsigned int index = ancestors.index_of(vertex);
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
			unsigned int index = ancestors.index_of(parent);
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
	unsigned int index = descendants.index_of(vertex);
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
			unsigned int index = descendants.index_of(child);
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

template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
struct set_info
{
	typedef typename ProofCalculus::Proof Proof;

	unsigned int set_size;
	Proof* size_axiom;
	hash_set<unsigned int> descendants;
	array<unsigned int> elements; /* NOTE: this does not contain the elements of the descendants of this set */

	static inline void free(set_info& info) {
		core::free(info.descendants);
		core::free(info.elements);
		core::free(*info.size_axiom);
		if (info.size_axiom->reference_count == 0)
			core::free(info.size_axiom);
		info.size_axiom = NULL;
	}

	inline Formula* set_formula() {
		return size_axiom->formula->binary.left->binary.right->quantifier.operand;
	}

	/* NOTE: this function does not check the consistency of the new size */
	inline void change_size(unsigned int new_size) {
#if !defined(NDEBUG)
		if (size_axiom->children.length != 0 && new_size != set_size)
			fprintf(stderr, "set_info.change_size WARNING: Attempted to change the size of a fixed set.\n");
#endif
		set_size = new_size;
		size_axiom->formula->binary.right->integer = new_size;
	}

	/* NOTE: this function does not check the consistency of the new size */
	inline bool set_size_axiom(Proof* axiom) {
		if (axiom == size_axiom) return true;
		if (size_axiom->children.length > 0) return false;
		core::free(*size_axiom);
		if (size_axiom->reference_count == 0)
			core::free(size_axiom);
		size_axiom = axiom;
		axiom->reference_count++;
		set_size = axiom->formula->binary.right->integer;
		return true;
	}
};

template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
inline bool init(
		set_info<BuiltInConstants, Formula, ProofCalculus>& info,
		unsigned int set_size, Formula* set_formula)
{
	info.set_size = set_size;
	Formula* size_axiom_formula = Formula::new_equals(Formula::new_atom(
			(unsigned int) BuiltInConstants::SIZE, Formula::new_lambda(1, set_formula)), Formula::new_int(set_size));
	if (size_axiom_formula == NULL) return false;
	set_formula->reference_count++;
	info.size_axiom = ProofCalculus::new_axiom(size_axiom_formula);
	free(*size_axiom_formula); if (size_axiom_formula->reference_count == 0) free(size_axiom_formula);
	if (info.size_axiom == NULL) return false;
	info.size_axiom->reference_count++;
	if (!hash_set_init(info.descendants, 16)) {
		free(*info.size_axiom); if (info.size_axiom->reference_count == 0) free(info.size_axiom);
		return false;
	} else if (!array_init(info.elements, 4)) {
		free(*info.size_axiom); if (info.size_axiom->reference_count == 0) free(info.size_axiom);
		free(info.descendants); return false;
	}
	return true;
}

template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
inline void on_free_set(unsigned int set_id,
		set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets)
{ }

template<typename Formula>
inline bool compute_new_set_size(
		const Formula* set_formula, unsigned int& out,
		unsigned int min_set_size, unsigned int max_set_size)
{
	out = (max_set_size == UINT_MAX) ? (min_set_size + 20) : ((min_set_size + max_set_size + 1) / 2);
	return true;
}

template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
struct set_reasoning
{
	typedef typename ProofCalculus::Proof Proof;

	struct changes {
		array_map<unsigned int, unsigned int> removed_sets;
	};

	extensional_set_graph<ProofCalculus> extensional_graph;
	intensional_set_graph intensional_graph;
	set_info<BuiltInConstants, Formula, ProofCalculus>* sets;

	unsigned int capacity;
	unsigned int set_count;

	hash_map<Formula, unsigned int> set_ids;
	hash_map<unsigned int, unsigned int> element_map;

	hash_multiset<unsigned int> symbols_in_formulas;

	set_reasoning() :
			extensional_graph(1024), intensional_graph(1024),
			capacity(1024), set_count(0), set_ids(2048),
			element_map(2048), symbols_in_formulas(256)
	{
		sets = (set_info<BuiltInConstants, Formula, ProofCalculus>*) malloc(sizeof(set_info<BuiltInConstants, Formula, ProofCalculus>) * capacity);
		if (sets == NULL) exit(EXIT_FAILURE);
		for (unsigned int i = 1; i < capacity; i++)
			sets[i].size_axiom = NULL;

		unsigned int empty_set_id;
		Formula* empty = Formula::new_false();
		if (empty == NULL
		 || !get_set_id(empty, empty_set_id))
			exit(EXIT_FAILURE);
		free(*empty); if (empty->reference_count == 0) free(empty);
		sets[empty_set_id].change_size(0);
	}

	~set_reasoning() {
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axiom != NULL) {
				extensional_graph.template free_set<false>(i);
				intensional_graph.free_set<false>(i);
				core::free(sets[i]);
			}
		} for (auto entry : set_ids)
			core::free(entry.key);
		core::free(sets);
	}

	bool ensure_capacity(unsigned int new_length) {
		if (new_length <= capacity) return true;
		unsigned int new_capacity = capacity;
		expand_capacity(new_capacity, new_length);

		if (!resize(sets, new_capacity) || !extensional_graph.resize(new_capacity) || !intensional_graph.resize(new_capacity))
			return false;
		for (unsigned int i = capacity; i < new_capacity; i++)
			sets[i].size_axiom = NULL; /* we use this field to mark which vertices are free */
		capacity = new_capacity;
		return true;
	}

	unsigned int get_next_free_set() const {
		for (unsigned int i = 1; i < capacity; i++)
			if (sets[i].size_axiom == NULL) return i;
		return capacity;
	}

	bool is_freeable(unsigned int set_id) const {
		return set_id > 1 && sets[set_id].size_axiom->reference_count == 1
			&& sets[set_id].elements.length == 0
			&& extensional_graph.vertices[set_id].children.size == 0
			&& extensional_graph.vertices[set_id].parents.size == 0;
	}

	inline void try_free_set(Formula* set_formula, unsigned int set_id) {
		if (is_freeable(set_id)) free_set(set_formula, set_id);
	}

	inline void try_free_set(Formula* set_formula) {
#if !defined(NDEBUG)
		bool contains;
		unsigned int set_id = set_ids.get(*set_formula, contains);
		if (!contains)
			fprintf(stderr, "set_reasoning.try_free_set WARNING: The given `set_formula` is not in `set_ids`.\n");
#else
		unsigned int set_id = set_ids.get(*set_formula);
#endif
		try_free_set(set_formula, set_id);
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
	bool new_set(Formula* set_formula, unsigned int& set_id, Args&&... visitor)
	{
		if (!new_set(set_id)) return false;

		/* initialize all intensional set relations */
		array<unsigned int> supersets(8), subsets(8);
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axiom == NULL) continue;
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
			if (!get_ancestors(intensional_graph, superset, ancestors)) {
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
		if (!init(sets[set_id], 1, set_formula)) {
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

		/* compute the descendants of `set_id` from its immediate children */
		for (unsigned int immediate_descendant : intensional_graph.vertices[set_id].children) {
			if (!sets[set_id].descendants.add_all(sets[immediate_descendant].descendants)) {
				free_set(set_id); return false;
			}
		}

		/* update the descendants set of all ancestors of `set_id` */
		hash_set<unsigned int> visited(32);
		array<unsigned int> stack(8);
		stack[stack.length++] = set_id;
		while (stack.length > 0) {
			unsigned int current = stack.pop();
			if (!sets[current].descendants.add(set_id)) {
				free_set(set_id); return false;
			}

			for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (visited.contains(parent)) continue;
				if (!visited.add(parent) || !stack.add(parent)) {
					free_set(set_id); return false;
				}
			} for (const auto& entry : extensional_graph.vertices[current].parents) {
				if (visited.contains(entry.key)) continue;
				if (!visited.add(entry.key) || !stack.add(entry.key)) {
					free_set(set_id); return false;
				}
			}
		}

		/* compute the upper bound and lower bound on the size of this new set */
		unsigned int min_set_size, max_set_size, initial_set_size;
		if (!set_size_bounds(set_id, min_set_size, max_set_size)
		 || !compute_new_set_size(set_formula, initial_set_size, min_set_size, max_set_size, std::forward<Args>(visitor)...))
		{
			free_set(set_id); return false;
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

	bool free_set(unsigned int set_id)
	{
#if !defined(NDEBUG)
		if (extensional_graph.vertices[set_id].parents.size > 0
		 || extensional_graph.vertices[set_id].children.size > 0)
			fprintf(stderr, "set_reasoning.free_set WARNING: This set has extensional edges.\n");
#endif

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

	inline bool free_set(Formula* formula, unsigned int set_id)
	{
		array_multiset<unsigned int> symbols(16);
		if (!get_constants(*formula, symbols)) return false;
		symbols_in_formulas.subtract<true>(symbols);

		bool contains;
		unsigned int bucket = set_ids.table.index_of(*formula, contains);
		core::free(set_ids.table.keys[bucket]);
		set_ids.remove_at(bucket);
		return free_set(set_id);
	}

	template<typename... Args>
	inline bool get_set_id(Formula* formula, unsigned int& set_id, bool& is_new, Args&&... visitor) {
		bool contains; unsigned int bucket;
		if (!set_ids.check_size()) return false;
		set_id = set_ids.get(*formula, contains, bucket);
		if (!contains) {
			array_multiset<unsigned int> symbols(16);
			if (!get_constants(*formula, symbols)) return false;
			symbols_in_formulas.add(symbols);

			if (!init(set_ids.table.keys[bucket], *formula)) {
				return false;
			} else if (!new_set(formula, set_id, std::forward<Args>(visitor)...)) {
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
	inline bool get_set_id(Formula* formula, unsigned int& set_id, Args&&... visitor) {
		bool is_new;
		return get_set_id(formula, set_id, is_new, std::forward<Args>(visitor)...);
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

	struct subgraph_formula_view {
		set_info<BuiltInConstants, Formula, ProofCalculus>* sets;
		unsigned int* indices;
		unsigned int length;

		subgraph_formula_view(
				set_info<BuiltInConstants, Formula, ProofCalculus>* sets,
				const hash_set<unsigned int>& vertices) : sets(sets), length(vertices.size)
		{
			indices = (unsigned int*) malloc(sizeof(unsigned int) * vertices.size);
			unsigned int index = 0;
			for (unsigned int vertex : vertices)
				indices[index++] = vertex;
		}

		~subgraph_formula_view() { free(indices); }

		Formula* operator[] (unsigned int index) const {
			return sets[indices[index]].set_formula();
		}

		inline unsigned int size() const {
			return length;
		}
	};

	bool uncontract_component(unsigned int contracted_set,
			const hash_set<unsigned int>& connected_component)
	{
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
					extensional_graph.vertices[entry.key].children.remove(contracted_set);
					extensional_graph.vertices[entry.key].children.put(member, extensional_graph.vertices[member].parents.get(entry.key));
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
					extensional_graph.vertices[entry.key].parents.remove(contracted_set);
					extensional_graph.vertices[entry.key].parents.put(member, extensional_graph.vertices[member].children.get(entry.key));
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
		free(sets[contracted_set]); if (contracted_set == set_count + 2) set_count--;
		return true;
	}

	bool contract_component(
			unsigned int set, unsigned int& contracted_set,
			const hash_set<unsigned int>& connected_component,
			bool& is_fixed)
	{
		is_fixed = false;

		Formula* conjunction = Formula::new_and(subgraph_formula_view(sets, connected_component));
		if (conjunction == NULL) return false;
		for (unsigned int member : connected_component)
			sets[member].set_formula()->reference_count++;

		if (!new_set(contracted_set)) {
			free(*conjunction); if (conjunction->reference_count == 0) free(conjunction);
			return false;
		} else if (!init(sets[contracted_set], sets[set].set_size, conjunction)) {
			free(*conjunction); if (conjunction->reference_count == 0) free(conjunction);
			intensional_graph.template free_set<true>(contracted_set);
			extensional_graph.template free_set<true>(contracted_set);
			return false;
		}
		free(*conjunction); if (conjunction->reference_count == 0) free(conjunction);
		if (contracted_set == set_count + 1) set_count++;

		if (!sets[contracted_set].descendants.add_all(sets[set].descendants)) {
			free(*conjunction); if (conjunction->reference_count == 0) free(conjunction);
			intensional_graph.template free_set<true>(contracted_set);
			extensional_graph.template free_set<true>(contracted_set);
			free(sets[contracted_set]); if (contracted_set == set_count + 2) set_count--;
			return false;
		}

		for (unsigned int member : connected_component) {
			if (!is_unfixed(member))
				is_fixed = true;
			for (const auto& entry : extensional_graph.vertices[member].parents) {
				if (!connected_component.contains(entry.key)) {
					extensional_graph.vertices[entry.key].children.remove(member);
					extensional_graph.vertices[entry.key].children.put(contracted_set, NULL);
					if (!extensional_graph.vertices[contracted_set].parents.put(entry.key, NULL)) {
						uncontract_component(contracted_set, connected_component);
						return false;
					}
				}
			} for (unsigned int parent : intensional_graph.vertices[member].parents) {
				if (!connected_component.contains(parent)) {
					unsigned int index = intensional_graph.vertices[parent].children.index_of(member);
					intensional_graph.vertices[parent].children.remove(index);
					if (!intensional_graph.vertices[parent].children.contains(contracted_set))
						intensional_graph.vertices[parent].children.add(contracted_set);
					if (!intensional_graph.vertices[contracted_set].parents.contains(parent)
					 && !intensional_graph.vertices[contracted_set].parents.add(parent)) {
						uncontract_component(contracted_set, connected_component);
						return false;
					}
				}
			} for (const auto& entry : extensional_graph.vertices[member].children) {
				if (!connected_component.contains(entry.key)) {
					extensional_graph.vertices[entry.key].parents.remove(member);
					extensional_graph.vertices[entry.key].parents.put(contracted_set, NULL);
					if (!extensional_graph.vertices[contracted_set].children.put(entry.key, NULL)) {
						uncontract_component(contracted_set, connected_component);
						return false;
					}
				}
			} for (unsigned int child : intensional_graph.vertices[member].children) {
				if (!connected_component.contains(child)) {
					unsigned int index = intensional_graph.vertices[child].parents.index_of(member);
					intensional_graph.vertices[child].parents.remove(index);
					if (!intensional_graph.vertices[child].parents.contains(contracted_set))
						intensional_graph.vertices[child].parents.add(contracted_set);
					if (!intensional_graph.vertices[contracted_set].children.contains(child)
					 && !intensional_graph.vertices[contracted_set].children.add(child)) {
						uncontract_component(contracted_set, connected_component);
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
					uncontract_component(contracted_set, connected_component);
					return false;
				}
			} for (const auto& entry : extensional_graph.vertices[ancestor].parents) {
				if (visited.contains(entry.key)) continue;
				if (!visited.add(entry.key) || !stack.add(entry.key)) {
					uncontract_component(contracted_set, connected_component);
					return false;
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
		if (!get_strongly_connected_component(set, connected_component)) {
			return false;
		} else if (connected_component.size > 1) {
			/* perform "surgery" on the graph by creating a new "virtual"
			   vertex that represents the strongly-connected component with all
			   of its internal edges contracted */
			if (!contract_component(set, set, connected_component, is_fixed))
				return false;
		} else {
			is_fixed = !is_unfixed(set);
		}

		if (!stack.ensure_capacity(stack.length + 1)) {
			if (connected_component.size > 1)
				uncontract_component(set, connected_component);
			return false;
		} if (stack.contains(set)) {
			/* we have previously requested to change the size of `set` */
			graph_changed = false; return true;
		}
		stack[stack.length++] = set;

		unsigned int* clique = NULL; unsigned int clique_count, ancestor_of_clique, upper_bound = 0;
		if (is_fixed) {
			upper_bound = sets[set].set_size;
		} else if (!get_size_upper_bound(set, upper_bound, clique, clique_count, ancestor_of_clique)) {
			if (connected_component.size > 1)
				uncontract_component(set, connected_component);
			return false;
		}
		if (requested_size <= upper_bound) {
			/* the bound is already satisfied */
			free(clique); stack.length--;
			graph_changed = true;
			if (connected_component.size > 1) {
				for (unsigned int member : connected_component)
					sets[member].change_size(requested_size);
				return uncontract_component(set, connected_component);
			} else {
				sets[set].change_size(requested_size);
				return true;
			}
		} else if (sets[set].set_size < upper_bound) {
			free(clique); stack.length--;
			graph_changed = true;
			if (connected_component.size > 1) {
				for (unsigned int member : connected_component)
					sets[member].change_size(upper_bound);
				return uncontract_component(set, connected_component);
			} else {
				sets[set].change_size(upper_bound);
				return true;
			}
		}

		if (!is_fixed) {
			bool child_graph_changed = false;
			if (!increase_set_size(ancestor_of_clique, sets[ancestor_of_clique].set_size + (requested_size - upper_bound), stack, child_graph_changed)) {
				if (connected_component.size > 1)
					uncontract_component(set, connected_component);
				free(clique); return false;
			} else if (child_graph_changed) {
				/* the graph has been changed */
				free(clique); stack.length--;
				graph_changed = true;
				if (connected_component.size > 1)
					return uncontract_component(set, connected_component);
				else return true;
			}

			for (unsigned int i = 0; i < clique_count; i++) {
				unsigned int requested_child_size;
				if (sets[clique[i]].set_size > requested_size - upper_bound)
					requested_child_size = sets[clique[i]].set_size - (requested_size - upper_bound);
				else requested_child_size = 0;
				if (!decrease_set_size(clique[i], requested_child_size, stack, child_graph_changed)) {
					if (connected_component.size > 1)
						uncontract_component(set, connected_component);
					free(clique); return false;
				} else if (child_graph_changed) {
					/* the graph has been changed */
					free(clique); stack.length--;
					graph_changed = true;
					if (connected_component.size > 1)
						uncontract_component(set, connected_component);
					else return true;
				}
			}
		}

		/* we were unable to change any bounds */
		free(clique); stack.length--;
		graph_changed = false;
		if (connected_component.size > 1)
			return uncontract_component(set, connected_component);
		else return true;
	}

	bool decrease_set_size(
			unsigned int set, unsigned int requested_size,
			array<unsigned int>& stack, bool& graph_changed)
	{
		/* first check if `set` is part of a strongly-connected component */
		bool is_fixed;
		hash_set<unsigned int> connected_component(16);
		if (!get_strongly_connected_component(set, connected_component)) {
			return false;
		} else if (connected_component.size > 1) {
			/* perform "surgery" on the graph by creating a new "virtual"
			   vertex that represents the strongly-connected component with all
			   of its internal edges contracted */
			if (!contract_component(set, set, connected_component, is_fixed))
				return false;
		} else {
			is_fixed = !is_unfixed(set);
		}

		if (!stack.ensure_capacity(stack.length + 1)) {
			if (connected_component.size > 1)
				uncontract_component(set, connected_component);
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
				uncontract_component(set, connected_component);
			return false;
		}
		if (requested_size >= lower_bound) {
			/* the bound is already satisfied */
			free(clique); stack.length--;
			graph_changed = true;
			if (connected_component.size > 1) {
				for (unsigned int member : connected_component)
					sets[member].change_size(requested_size);
				return uncontract_component(set, connected_component);
			} else {
				sets[set].change_size(requested_size);
				return true;
			}
		} else if (sets[set].set_size > lower_bound) {
			free(clique); stack.length--;
			graph_changed = true;
			if (connected_component.size > 1) {
				for (unsigned int member : connected_component)
					sets[member].change_size(lower_bound);
				return uncontract_component(set, connected_component);
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
							uncontract_component(set, connected_component);
						free(clique); return false;
					} else if (child_graph_changed) {
						/* the graph has been changed */
						free(clique); stack.length--;
						graph_changed = true;
						if (connected_component.size > 1)
							return uncontract_component(set, connected_component);
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
						uncontract_component(set, connected_component);
					free(clique); return false;
				} else if (child_graph_changed) {
					/* the graph has been changed */
					free(clique); stack.length--;
					graph_changed = true;
					if (connected_component.size > 1)
						return uncontract_component(set, connected_component);
					else return true;
				}

				for (unsigned int i = 0; i < clique_count; i++) {
					unsigned int requested_child_size;
					if (sets[clique[i]].set_size > requested_ancestor_size - sets[ancestor_of_clique].set_size)
						requested_child_size = sets[clique[i]].set_size - (requested_ancestor_size - sets[ancestor_of_clique].set_size);
					else requested_child_size = 0;
					if (!decrease_set_size(clique[i], requested_child_size, stack, child_graph_changed)) {
						if (connected_component.size > 1)
							uncontract_component(set, connected_component);
						free(clique); return false;
					} else if (child_graph_changed) {
						/* the graph has been changed */
						free(clique); stack.length--;
						graph_changed = true;
						if (connected_component.size > 1)
							return uncontract_component(set, connected_component);
						else return true;
					}
				}
			}
		}

		/* we were unable to change any bounds */
		free(clique); stack.length--;
		graph_changed = false;
		if (connected_component.size > 1)
			return uncontract_component(set, connected_component);
		else return true;
	}

	/* NOTE: this function does not check for consistency */
	bool add_subset_axiom(Proof* axiom)
	{
		Formula* antecedent = axiom->formula->quantifier.operand->binary.left;
		Formula* consequent = axiom->formula->quantifier.operand->binary.right;

		if (!set_ids.check_size(set_ids.table.size + 2)) return false;

		unsigned int antecedent_set, consequent_set;
		bool is_antecedent_new, is_consequent_new;
		if (!get_set_id(antecedent, antecedent_set, is_antecedent_new)) {
			return false;
		} else if (!get_set_id(consequent, consequent_set, is_consequent_new)) {
			if (is_freeable(antecedent_set)) free_set(antecedent, antecedent_set);
			return false;
		}

#if !defined(NDEBUG)
		if (consequent_set == antecedent_set)
			fprintf(stderr, "set_reasoning.add_subset_axiom WARNING: `consequent` and `antecedent` are the same set.\n");
#endif

		if (!extensional_graph.add_edge(consequent_set, antecedent_set, axiom)) {
			/* if either the antecedent or consequent sets have no references, free them */
			if (is_freeable(consequent_set)) free_set(consequent, consequent_set);
			if (is_freeable(antecedent_set)) free_set(antecedent, antecedent_set);
			return false;
		}
		return true;
	}

	template<bool ResolveInconsistencies, typename... Args>
	Proof* get_subset_axiom(Formula* antecedent, Formula* consequent, Args&&... visitor)
	{
		if (!set_ids.check_size(set_ids.table.size + 2)) return NULL;

		unsigned int antecedent_set, consequent_set;
		bool is_antecedent_new, is_consequent_new;
		if (!get_set_id(antecedent, antecedent_set, is_antecedent_new, std::forward<Args>(visitor)...)) {
			return NULL;
		} else if (!get_set_id(consequent, consequent_set, is_consequent_new, std::forward<Args>(visitor)...)) {
			if (is_freeable(antecedent_set)) free_set(antecedent, antecedent_set);
			return NULL;
		}

#if !defined(NDEBUG)
		if (consequent_set == antecedent_set)
			fprintf(stderr, "set_reasoning.get_subset_axiom WARNING: `consequent` and `antecedent` are the same set.\n");
#endif

		bool new_edge;
		Proof* axiom = extensional_graph.get_edge(consequent_set, antecedent_set, consequent, antecedent, new_edge);
		if (axiom == NULL) {
			/* if either the antecedent or consequent sets have no references, free them */
			if (is_freeable(consequent_set)) free_set(consequent, consequent_set);
			if (is_freeable(antecedent_set)) free_set(antecedent, antecedent_set);
			return NULL;
		}

		if (!new_edge) return axiom;

		/* check that the new edge does not create any inconsistencies */
		for (unsigned int descendant : sets[antecedent_set].descendants) {
			if (sets[descendant].set_size > 0 && are_disjoint(consequent_set, descendant)) {
				if (ResolveInconsistencies) {
					while (sets[descendant].set_size > 0) {
						bool graph_changed;
						array<unsigned int> stack(8);
						if (!decrease_set_size(descendant, 0, stack, graph_changed)) {
							extensional_graph.remove_edge(consequent_set, antecedent_set);
							if (is_freeable(consequent_set)) free_set(consequent, consequent_set);
							if (is_freeable(antecedent_set)) free_set(antecedent, antecedent_set);
							return NULL;
						} else if (graph_changed) continue;

						/* we were unable to change the graph, so it cannot be made consistent (as far as we can tell) */
						extensional_graph.remove_edge(consequent_set, antecedent_set);
						if (is_freeable(consequent_set)) free_set(consequent, consequent_set);
						if (is_freeable(antecedent_set)) free_set(antecedent, antecedent_set);
						return NULL;
					}
				} else {
					extensional_graph.remove_edge(consequent_set, antecedent_set);
					if (is_freeable(consequent_set)) free_set(consequent, consequent_set);
					if (is_freeable(antecedent_set)) free_set(antecedent, antecedent_set);
					return NULL;
				}
			}
		}
		while (true) {
			unsigned int* clique = NULL; unsigned int clique_count; unsigned int ancestor_of_clique;
			if (sets[antecedent_set].set_size == 0) {
				break;
			} else if (sets[consequent_set].set_size == 0) {
				bool graph_changed;
				array<unsigned int> stack(8);
				if (!decrease_set_size(antecedent_set, 0, stack, graph_changed)) {
					extensional_graph.remove_edge(consequent_set, antecedent_set);
					if (is_freeable(consequent_set)) free_set(consequent, consequent_set);
					if (is_freeable(antecedent_set)) free_set(antecedent, antecedent_set);
					return NULL;
				} else if (graph_changed) continue;

				/* we were unable to change the graph, so it cannot be made consistent (as far as we can tell) */
				extensional_graph.remove_edge(consequent_set, antecedent_set);
				if (is_freeable(consequent_set)) free_set(consequent, consequent_set);
				if (is_freeable(antecedent_set)) free_set(antecedent, antecedent_set);
				return NULL;
			} else if (!find_largest_disjoint_clique_with_set(*this, antecedent_set, consequent_set, clique, clique_count, ancestor_of_clique, 1)) {
				extensional_graph.remove_edge(consequent_set, antecedent_set);
				if (is_freeable(consequent_set)) free_set(consequent, consequent_set);
				if (is_freeable(antecedent_set)) free_set(antecedent, antecedent_set);
			}

			if (clique == NULL) break;

			if (ResolveInconsistencies) {
				/* try increasing the size of the ancestor */
				unsigned int clique_size = 0;
				for (unsigned int i = 0; i < clique_count; i++)
					clique_size += sets[clique[i]].set_size;

				bool graph_changed;
				array<unsigned int> stack(8);
				if (!increase_set_size(ancestor_of_clique, clique_size, stack, graph_changed)) {
					extensional_graph.remove_edge(consequent_set, antecedent_set);
					if (is_freeable(consequent_set)) free_set(consequent, consequent_set);
					if (is_freeable(antecedent_set)) free_set(antecedent, antecedent_set);
					free(clique); return NULL;
				} else if (graph_changed) { free(clique); continue; }

				for (unsigned int i = 0; i < clique_count; i++) {
					unsigned int requested_size;
					if (sets[ancestor_of_clique].set_size > clique_size - sets[clique[i]].set_size)
						requested_size = sets[ancestor_of_clique].set_size - (clique_size - sets[clique[i]].set_size);
					else requested_size = 0;

					if (!decrease_set_size(clique[i], requested_size, stack, graph_changed)) {
						extensional_graph.remove_edge(consequent_set, antecedent_set);
						if (is_freeable(consequent_set)) free_set(consequent, consequent_set);
						if (is_freeable(antecedent_set)) free_set(antecedent, antecedent_set);
						free(clique); return NULL;
					} else if (graph_changed) break;
				}

				if (graph_changed) { free(clique); continue; }

				/* we were unable to change the graph, so it cannot be made consistent (as far as we can tell) */
				extensional_graph.remove_edge(consequent_set, antecedent_set);
				if (is_freeable(consequent_set)) free_set(consequent, consequent_set);
				if (is_freeable(antecedent_set)) free_set(antecedent, antecedent_set);
				free(clique); return NULL;
			} else {
				extensional_graph.remove_edge(consequent_set, antecedent_set);
				if (is_freeable(consequent_set)) free_set(consequent, consequent_set);
				if (is_freeable(antecedent_set)) free_set(antecedent, antecedent_set);
				free(clique); return NULL;
			}
		}

		/* update the descendants of `consequent_set` and all its ancestors */
		array<pair<unsigned int, unsigned int>> stack(8);
		stack[stack.length++] = {consequent_set, antecedent_set};
		while (stack.length > 0) {
			pair<unsigned int, unsigned int> current = stack.pop();
			unsigned int old_size = sets[current.key].descendants.size;
			if (!sets[current.key].descendants.add_all(sets[current.value].descendants)) {
				remove_subset_relation(antecedent_set, consequent_set, antecedent, consequent);
				return NULL;
			}
			if (old_size == sets[current.key].descendants.size) continue;

			for (unsigned int parent : intensional_graph.vertices[current.key].parents) {
				if (!stack.add({parent, current.key})) {
					remove_subset_relation(antecedent_set, consequent_set, antecedent, consequent);
					return NULL;
				}
			} for (const auto& entry : extensional_graph.vertices[current.key].parents) {
				if (!stack.add({entry.key, current.key})) {
					remove_subset_relation(antecedent_set, consequent_set, antecedent, consequent);
					return NULL;
				}
			}
		}
		return axiom;
	}

	bool remove_subset_relation(
			unsigned int antecedent_set, unsigned int consequent_set,
			Formula* antecedent, Formula* consequent)
	{
#if !defined(NDEBUG)
		if (consequent_set == antecedent_set)
			fprintf(stderr, "set_reasoning.remove_subset_relation WARNING: `consequent` and `antecedent` are the same set.\n");
#endif

		extensional_graph.remove_edge(consequent_set, antecedent_set);

		array<unsigned int> stack(8);
		stack[stack.length++] = consequent_set;
		while (stack.length > 0) {
			unsigned int current = stack.pop();
			unsigned int old_size = sets[current].descendants.size;
			sets[current].descendants.clear();
			sets[current].descendants.add(current);
			for (unsigned int child : intensional_graph.vertices[current].children)
				if (!sets[current].descendants.add_all(sets[child].descendants)) return false;
			for (const auto& entry : extensional_graph.vertices[current].children)
				if (!sets[current].descendants.add_all(sets[entry.key].descendants)) return false;
			if (old_size == sets[current].descendants.size) continue;

			for (const auto& entry : extensional_graph.vertices[current].parents)
				if (!stack.add(entry.key)) return false;
			for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (extensional_graph.vertices[current].parents.contains(parent)) continue;
				if (!stack.add(parent)) return false;
			}
		}

		/* if either the antecedent or consequent sets have no references, free them */
		if (is_freeable(consequent_set) && !free_set(consequent, consequent_set)) return false;
		if (is_freeable(antecedent_set) && !free_set(antecedent, antecedent_set)) return false;
		return true;
	}

	bool free_subset_axiom(Formula* antecedent, Formula* consequent)
	{
		antecedent->reference_count++;
		consequent->reference_count++;

#if !defined(NDEBUG)
		bool contains;
		unsigned int antecedent_set = set_ids.get(*antecedent, contains);
		if (!contains)
			fprintf(stderr, "set_reasoning.try_free_subset_axiom WARNING: No such set for given antecedent.\n");
		unsigned int consequent_set = set_ids.get(*consequent, contains);
		if (!contains)
			fprintf(stderr, "set_reasoning.try_free_subset_axiom WARNING: No such set for given consequent.\n");
#else
		unsigned int antecedent_set = set_ids.get(*antecedent);
		unsigned int consequent_set = set_ids.get(*consequent);
#endif

		bool success = remove_subset_relation(antecedent_set, consequent_set, antecedent, consequent);
		free(*antecedent); if (antecedent->reference_count == 0) free(antecedent);
		free(*consequent); if (consequent->reference_count == 0) free(consequent);
		return success;
	}

	inline bool free_subset_axiom(Proof* subset_axiom)
	{
		Formula* antecedent = subset_axiom->formula->quantifier.operand->binary.left;
		Formula* consequent = subset_axiom->formula->quantifier.operand->binary.right;
		return free_subset_axiom(antecedent, consequent);
	}

	template<bool ResolveInconsistencies>
	Proof* get_size_axiom(unsigned int set_id, unsigned int new_size)
	{
		bool is_fixed = !is_unfixed(set_id);
		if (new_size > sets[set_id].set_size) {
			unsigned int upper_bound = 0;
			if (is_fixed) {
				upper_bound = sets[set_id].set_size;
			} else if (!get_size_upper_bound(set_id, upper_bound))
				return NULL;
			if (new_size <= upper_bound) {
				sets[set_id].change_size(new_size);
				return sets[set_id].size_axiom;
			}

			if (!ResolveInconsistencies || is_fixed) return NULL;

			while (true) {
				/* compute the upper bound on the size of this set; if the new size
				   violates this bound, change the sizes of other sets to increase the bound */
				array<unsigned int> stack(8); bool graph_changed;
				if (!increase_set_size(set_id, new_size, stack, graph_changed))
					return NULL;
				if (sets[set_id].set_size == new_size)
					break;
				if (!graph_changed) return NULL;
			}
		} else if (new_size < sets[set_id].set_size) {
			unsigned int lower_bound = 0;
			if (is_fixed) {
				lower_bound = sets[set_id].set_size;
			} else if (!get_size_lower_bound(set_id, lower_bound))
				return NULL;
			if (new_size >= lower_bound) {
				sets[set_id].change_size(new_size);
				return sets[set_id].size_axiom;
			}

			if (!ResolveInconsistencies || is_fixed) return NULL;

			while (true) {
				/* compute the lower bound on the size of this set; if the new size
				   violates this bound, change the sizes of other sets to increase the bound */
				array<unsigned int> stack(8); bool graph_changed;
				if (!decrease_set_size(set_id, new_size, stack, graph_changed))
					return NULL;
				if (sets[set_id].set_size == new_size)
					break;
				if (!graph_changed) return NULL;
			}
		}
		return sets[set_id].size_axiom;
	}

	template<bool ResolveInconsistencies>
	inline Proof* get_size_axiom(Formula* formula, unsigned int new_size) {
		unsigned int set_id;
		if (!get_set_id(formula, set_id))
			return NULL;
		return get_size_axiom<ResolveInconsistencies>(set_id, new_size);
	}

	inline bool get_provable_elements(unsigned int set_id,
			hash_set<unsigned int>& provable_elements) const
	{
		for (unsigned int descendant : sets[set_id].descendants) {
			for (unsigned int element : sets[descendant].elements)
				if (!provable_elements.add(element)) return false;
		}
		return true;
	}

	bool get_size_lower_bound(unsigned int set_id, unsigned int& lower_bound,
			unsigned int*& clique, unsigned int& clique_count,
			unsigned int& ancestor_of_clique) const
	{
#if !defined(NDEBUG)
		if (!is_unfixed(set_id))
			fprintf(stderr, "get_size_lower_bound WARNING: Set with ID %u is fixed. This function should not be called on fixed sets.\n", set_id);
#endif

		/* first compute the number of elements that provably belong to this set */
		hash_set<unsigned int> provable_elements(16);
		if (!get_provable_elements(set_id, provable_elements))
			return false;
		lower_bound = provable_elements.size;

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
			bool can_be_zero = true;
			for (unsigned int i = 1; can_be_zero && i < set_count + 1; i++) {
				if (sets[i].size_axiom == NULL) continue;
				for (unsigned int j = i + 1; j < set_count + 1; j++) {
					if (sets[j].size_axiom == NULL) continue;
					if (are_newly_disjoint(sets[i].set_formula(), sets[j].set_formula(), set_id)) {
						/* check if any ancestor of `i` and `j` has a set size smaller than its lower bound */
						if (!find_largest_disjoint_clique_with_edge(*this, i, j, clique, clique_count, ancestor_of_clique, INT_MIN, 0)) {
							sets[set_id].set_size = old_size;
							return false;
						}
						if (clique != NULL) {
							unsigned int clique_size = 0;
							for (unsigned int i = 0; i < clique_count; i++)
								clique_size += sets[clique[i]].set_size;
							if (clique_size > sets[ancestor_of_clique].set_size) {
								lower_bound = 1;
								can_be_zero = false;
								break;
							} else {
								free(clique);
								clique = NULL;
								clique_count = 0;
								ancestor_of_clique = 0;
							}
						}
					}
				}
			}
			sets[set_id].set_size = old_size;
		} else {
			ancestor_of_clique = 0;
		}
		return true;
	}

	inline bool get_size_lower_bound(unsigned int set_id, unsigned int& lower_bound) const
	{
		unsigned int* clique = NULL; unsigned int clique_count, ancestor_of_clique;
		if (!get_size_lower_bound(set_id, lower_bound, clique, clique_count, ancestor_of_clique))
			return false;
		if (clique != NULL) free(clique);
		return true;
	}

	bool get_size_upper_bound(unsigned int set_id, unsigned int& upper_bound,
			unsigned int*& clique, unsigned int& clique_count, unsigned int& ancestor_of_clique) const
	{
#if !defined(NDEBUG)
		if (!is_unfixed(set_id))
			fprintf(stderr, "get_size_upper_bound WARNING: Set with ID %u is fixed. This function should not be called on fixed sets.\n", set_id);
#endif

		if (!find_largest_disjoint_clique_with_set(*this, set_id, clique, clique_count, ancestor_of_clique, INT_MIN))
			return false;
		if (clique == NULL) {
			upper_bound = UINT_MAX; /* the upper bound is infinite */
			return true;
		}

		unsigned int index = index_of(set_id, clique, clique_count);
		clique[index] = clique[--clique_count];

		upper_bound = sets[ancestor_of_clique].set_size;
		for (unsigned int i = 0; i < clique_count; i++)
			upper_bound -= sets[clique[i]].set_size;
		return true;
	}

	inline bool get_size_upper_bound(unsigned int set_id, unsigned int& upper_bound) const
	{
		unsigned int* clique = NULL; unsigned int clique_count; unsigned int ancestor_of_clique;
		if (!get_size_upper_bound(set_id, upper_bound, clique, clique_count, ancestor_of_clique))
			return false;
		if (clique != NULL) free(clique);
		return true;
	}

	template<bool ResolveInconsistencies, typename... Args>
	bool move_element_to_set(unsigned int element, Formula* new_set_formula, Args&&... visitor)
	{
		if (!element_map.check_size()) return false;
		bool contains; unsigned int bucket;
		unsigned int old_set_id = element_map.get(element, contains, bucket);
		if (!contains) old_set_id = 0;

#if !defined(NDEBUG)
		typedef typename Formula::Type FormulaType;

		if (new_set_formula->type == FormulaType::FALSE) {
			fprintf(stderr, "move_element_to_set WARNING: `new_set_formula` is false.\n");
			return false;
		}
#endif

		unsigned int new_set_id;
		if (!get_set_id(new_set_formula, new_set_id, std::forward<Args>(visitor)...))
			return false;
		else if (new_set_id == old_set_id)
			return true;

		if (old_set_id != 0) {
			unsigned int index = sets[old_set_id].elements.index_of(element);
#if !defined(NDEBUG)
			if (index == sets[old_set_id].elements.length)
				fprintf(stderr, "set_reasoning.move_element_to_set WARNING: The old set does not contain the given element.\n");
			else sets[old_set_id].elements.remove(index);
#else
			sets[old_set_id].elements.remove(index);
#endif
		}

		if (!sets[new_set_id].elements.ensure_capacity(sets[new_set_id].elements.length + 1)) {
			if (is_freeable(new_set_id)) free_set(new_set_formula, new_set_id);
			return false;
		}
#if !defined(NDEBUG)
		unsigned int index = sets[new_set_id].elements.index_of(element);
		if (index < sets[new_set_id].elements.length)
			fprintf(stderr, "set_reasoning.move_element_to_set WARNING: The new set already contains the given element.\n");
#endif
		sets[new_set_id].elements[sets[new_set_id].elements.length++] = element;

		/* retrieve the ancestors of `old_set_id` */
		array<unsigned int> stack(8);
		hash_set<unsigned int> old_ancestors(16);
		if (old_set_id != 0) {
			stack[stack.length++] = old_set_id;
			old_ancestors.add(old_set_id);
			while (stack.length > 0) {
				unsigned int current = stack.pop();
				for (unsigned int parent : intensional_graph.vertices[current].parents) {
					if (old_ancestors.contains(parent)) continue;
					if (!old_ancestors.add(parent) || !stack.add(parent)) {
						if (old_set_id != 0) sets[old_set_id].elements.add(element);
						sets[new_set_id].elements.length--;
						if (is_freeable(new_set_id)) free_set(new_set_formula, new_set_id);
						return false;
					}
				} for (const auto& entry : extensional_graph.vertices[current].parents) {
					if (old_ancestors.contains(entry.key)) continue;
					if (!old_ancestors.add(entry.key) || !stack.add(entry.key)) {
						if (old_set_id != 0) sets[old_set_id].elements.add(element);
						sets[new_set_id].elements.length--;
						if (is_freeable(new_set_id)) free_set(new_set_formula, new_set_id);
						return false;
					}
				}
			}
		}

		/* retrieve the new ancestors of `new_set_id` */
		hash_set<unsigned int> new_ancestors(16);
		stack[stack.length++] = new_set_id;
		new_ancestors.add(new_set_id);
		while (stack.length > 0) {
			unsigned int current = stack.pop();
			for (unsigned int parent : intensional_graph.vertices[current].parents) {
				if (old_ancestors.contains(parent) || new_ancestors.contains(parent)) continue;
				if (!new_ancestors.add(parent) || !stack.add(parent)) {
					if (old_set_id != 0) sets[old_set_id].elements.add(element);
					sets[new_set_id].elements.length--;
					if (is_freeable(new_set_id)) free_set(new_set_formula, new_set_id);
					return false;
				}
			} for (const auto& entry : extensional_graph.vertices[current].parents) {
				if (old_ancestors.contains(entry.key) || new_ancestors.contains(entry.key)) continue;
				if (!new_ancestors.add(entry.key) || !stack.add(entry.key)) {
					if (old_set_id != 0) sets[old_set_id].elements.add(element);
					sets[new_set_id].elements.length--;
					if (is_freeable(new_set_id)) free_set(new_set_formula, new_set_id);
					return false;
				}
			}
		}

		/* check the new ancestors of `new_set_id` to make sure the size of each
		   set is at least as large as the number of provable elements of each set */
		for (unsigned int new_ancestor : new_ancestors) {
			hash_set<unsigned int> provable_elements(16);
			if (!get_provable_elements(new_ancestor, provable_elements)) {
				if (old_set_id != 0) sets[old_set_id].elements.add(element);
				sets[new_set_id].elements.length--;
				if (is_freeable(new_set_id)) free_set(new_set_formula, new_set_id);
				return false;
			}

			if (sets[new_ancestor].set_size >= provable_elements.size) continue;

			while (true) {
				/* compute the upper bound on the size of this set; if the new size
					violates this bound, change the sizes of other sets to increase the bound */
				array<unsigned int> stack(8); bool graph_changed;
				if (!ResolveInconsistencies || !increase_set_size(new_ancestor, provable_elements.size, stack, graph_changed)) {
					if (old_set_id != 0) sets[old_set_id].elements.add(element);
					sets[new_set_id].elements.length--;
					if (is_freeable(new_set_id)) free_set(new_set_formula, new_set_id);
					return false;
				}
				if (sets[new_ancestor].set_size >= provable_elements.size)
					break;
				if (!graph_changed) {
					if (old_set_id != 0) sets[old_set_id].elements.add(element);
					sets[new_set_id].elements.length--;
					if (is_freeable(new_set_id)) free_set(new_set_formula, new_set_id);
					return false;
				}
			}
		}

		/* the set `old_set_id` might be freeable now */
		if (old_set_id == 0) {
			element_map.table.keys[bucket] = element;
			element_map.table.size++;
		}
		element_map.values[bucket] = new_set_id;
		if (is_freeable(old_set_id)) {
			on_free_set(old_set_id, *this, std::forward<Args>(visitor)...);
			free_set(sets[old_set_id].set_formula(), old_set_id);
		}
		return true;
	}

	template<typename... Args>
	bool move_element_to_superset(unsigned int element, Formula* new_set_formula, Args&&... visitor)
	{
		bool contains; unsigned int bucket;
		unsigned int old_set_id = element_map.get(element, contains, bucket);
#if !defined(NDEBUG)
		if (!contains) {
			fprintf(stderr, "set_reasoning.move_element_to_superset ERROR: The given element is not in any set.\n");
			return false;
		}
#endif

		unsigned int new_set_id;
		if (!get_set_id(new_set_formula, new_set_id, std::forward<Args>(visitor)...))
			return false;
		else if (new_set_id == old_set_id)
			return true;

#if !defined(NDEBUG)
		if (!sets[new_set_id].descendants.contains(old_set_id))
			fprintf(stderr, "set_reasoning.move_element_to_superset WARNING: The new set is not a superset of the old set.\n");
#endif

		unsigned int index = sets[old_set_id].elements.index_of(element);
#if !defined(NDEBUG)
		if (index == sets[old_set_id].elements.length)
			fprintf(stderr, "set_reasoning.move_element_to_superset WARNING: The old set does not contain the given element.\n");
		else sets[old_set_id].elements.remove(index);
#else
		sets[old_set_id].elements.remove(index);
#endif

		if (!sets[new_set_id].elements.ensure_capacity(sets[new_set_id].elements.length + 1)) {
			if (is_freeable(new_set_id)) free_set(new_set_formula, new_set_id);
			return false;
		}
#if !defined(NDEBUG)
		index = sets[new_set_id].elements.index_of(element);
		if (index < sets[new_set_id].elements.length)
			fprintf(stderr, "set_reasoning.move_element_to_superset WARNING: The new set already contains the given element.\n");
#endif
		sets[new_set_id].elements[sets[new_set_id].elements.length++] = element;

		/* `old_set_id` may now be freeable */
		element_map.values[bucket] = new_set_id;
		if (is_freeable(old_set_id)) {
			on_free_set(old_set_id, *this, std::forward<Args>(visitor)...);
			free_set(sets[old_set_id].set_formula(), old_set_id);
		}
		return true;
	}

	template<typename... Args>
	void remove_element(unsigned int element, Args&&... visitor)
	{
		bool contains; unsigned int bucket;
		unsigned int set_id = element_map.get(element, contains, bucket);
#if !defined(NDEBUG)
		if (!contains) {
			fprintf(stderr, "set_reasoning.remove_element ERROR: The given element is not in any set.\n");
			return;
		}
#endif

		unsigned int index = sets[set_id].elements.index_of(element);
#if !defined(NDEBUG)
		if (index == sets[set_id].elements.length) {
			fprintf(stderr, "set_reasoning.remove_element ERROR: The set does not contain the given element.\n");
			return;
		}
#endif
		sets[set_id].elements.remove(index);

		/* `set_id` may now be freeable */
		element_map.remove_at(bucket);
		if (is_freeable(set_id)) {
			on_free_set(set_id, *this, std::forward<Args>(visitor)...);
			free_set(sets[set_id].set_formula(), set_id);
		}
	}

	inline bool are_disjoint(Formula* first, Formula* second) const
	{
		Formula* intersection = intersect(first, second);

		array<unsigned int> stack(8);
		hash_set<unsigned int> visited(16);
		stack[stack.length++] = 1;
		visited.add(1);
		while (stack.length > 0) {
			unsigned int current = stack.pop();
			if (is_subset(intersection, sets[current].set_formula())) {
				free(*intersection); if (intersection->reference_count == 0) free(intersection);
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
		free(*intersection); if (intersection->reference_count == 0) free(intersection);
		return false;
	}

	inline bool are_disjoint(unsigned int first_set, unsigned int second_set) const
	{
		return are_disjoint(sets[first_set].set_formula(), sets[second_set].set_formula());
	}

	/**
	 * Returns `true` when `first` and `second` are disjoint iff the set with
	 * ID `set_id` is empty (and all other sets stay the same).
	 */
	inline bool are_newly_disjoint(Formula* first, Formula* second, unsigned int set_id) const
	{
		Formula* intersection = intersect(first, second);

		/* first check if the intersection belongs to `set_id` */
		if (!is_subset(intersection, sets[set_id].set_formula())) {
			free(*intersection); if (intersection->reference_count == 0) free(intersection);
			return false;
		}

		/* then check if the intersection belongs to any child */
		for (const auto& entry : extensional_graph.vertices[set_id].children) {
			if (is_subset(intersection, sets[entry.key].set_formula())) {
				free(*intersection); if (intersection->reference_count == 0) free(intersection);
				return false;
			}
		} for (unsigned int child : intensional_graph.vertices[set_id].children) {
			if (is_subset(intersection, sets[child].set_formula())) {
				free(*intersection); if (intersection->reference_count == 0) free(intersection);
				return false;
			}
		}

		free(*intersection); if (intersection->reference_count == 0) free(intersection);
		return true;
	}

	inline bool is_unfixed(Proof* size_axiom) const {
		return (size_axiom->children.length == 0);
	}

	inline bool is_unfixed(Proof* size_axiom, const hash_set<Proof*>& observations) const {
		return (size_axiom->children.length == 0
			&& !observations.contains(size_axiom));
	}

	inline bool is_unfixed(unsigned int set_id) const {
		return is_unfixed(sets[set_id].size_axiom);
	}

	inline bool is_unfixed(unsigned int set_id, const hash_set<Proof*>& observations) const {
		return is_unfixed(sets[set_id].size_axiom, observations);
	}

	bool get_unfixed_sets(array<unsigned int>& unfixed_set_ids, const hash_set<Proof*>& observations) {
		for (unsigned int i = 2; i < set_count + 1; i++) {
			if (sets[i].size_axiom == NULL || !is_unfixed(i, observations))
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
			if (sets[i].size_axiom != NULL) {
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

	template<typename Stream, typename... Printer>
	bool print_axioms(Stream& out, Printer&&... printer) const {
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axiom == NULL) continue;
			if (!print(*sets[i].size_axiom->formula, out, std::forward<Printer>(printer)...) || !print('\n', out))
				return false;
			for (const auto& entry : extensional_graph.vertices[i].children) {
				if (!print(*entry.value->formula, out, std::forward<Printer>(printer)...) || !print('\n', out))
					return false;
			}
		}
		return true;
	}

	bool are_elements_unique() const {
		hash_map<unsigned int, array<unsigned int>> all_elements(1024);
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axiom == NULL) continue;
			if (!all_elements.check_size(all_elements.table.size + sets[i].elements.length)) {
				for (auto entry : all_elements) free(entry.value);
				return false;
			}
			for (unsigned int element : sets[i].elements) {
				bool contains; unsigned int bucket;
				array<unsigned int>& containing_sets = all_elements.get(element, contains, bucket);
				if (!contains) {
					if (!array_init(containing_sets, 4)) {
						for (auto entry : all_elements) free(entry.value);
						return false;
					}
					all_elements.table.keys[bucket] = element;
					all_elements.table.size++;
				} if (!containing_sets.add(i)) {
					for (auto entry : all_elements) free(entry.value);
					return false;
				}
			}
		}

		bool success = true;
		for (const auto& entry : all_elements) {
			if (entry.value.length > 1) {
				print("set_reasoning.are_elements_unique WARNING: Detected element ", stderr);
				print(entry.key, stderr); print(" in multiple sets: ", stderr);
				print<unsigned int, empty_string, empty_string>(entry.value, stderr); print('\n', stderr);
				success = false;
			}
		}

		for (auto entry : all_elements) free(entry.value);
		return success;
	}

	void check_freeable_sets() const {
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axiom == NULL) continue;
			if (is_freeable(i))
				fprintf(stderr, "set_reasoning.check_freeable_sets: Set with ID %u is freeable.\n", i);
		}
	}

	bool are_set_sizes_valid() const {
		bool success = true;
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axiom == NULL) continue;
			unsigned int min_set_size, max_set_size;
			if (!set_size_bounds(i, min_set_size, max_set_size))
				return false;
			if (sets[i].set_size < min_set_size || sets[i].set_size > max_set_size) {
				fprintf(stderr, "set_reasoning.are_set_sizes_valid WARNING: Set with ID %u has size outside the bounds computed by `set_reasoning.set_size_bounds`.\n", i);
				success = false;
			}
		}
		return success;
	}

	void check_set_ids() const {
		for (unsigned int i = 1; i < set_count + 1; i++) {
			if (sets[i].size_axiom == NULL) continue;

			bool contains;
			unsigned int set_id = set_ids.get(*sets[i].set_formula(), contains);
			if (!contains) {
				fprintf(stderr, "set_reasoning.check_set_ids WARNING: Set with ID %u has set formula '", i);
				print(*sets[i].set_formula(), stderr); fprintf(stderr, "' that does not exist in the `set_ids` map.\n");
			} else if (set_id != i) {
				fprintf(stderr, "set_reasoning.check_set_ids WARNING: Set with ID %u has set formula '", i);
				print(*sets[i].set_formula(), stderr); fprintf(stderr, "' that is mapped to %u in the `set_ids` map.\n", set_id);
			}
		}

		for (const auto& entry : set_ids) {
			if (sets[entry.value].size_axiom == NULL) {
				fprintf(stderr, "set_reasoning.check_set_ids WARNING: `set_ids` map contains an entry from formula '");
				print(entry.key, stderr);
				fprintf(stderr, "' to %u, but the set with ID %u doesn't exist (the size axiom is null).\n", entry.value, entry.value);
			} else if (*sets[entry.value].set_formula() != entry.key) {
				fprintf(stderr, "set_reasoning.check_set_ids WARNING: `set_ids` map contains an entry from formula '");
				print(entry.key, stderr); fprintf(stderr, "' to %u, but the set with ID %u has set formula '", entry.value, entry.value);
				print(*sets[entry.value].set_formula(), stderr); fprintf(stderr, "'.\n");
			}
		}
	}
};

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

template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
inline bool init(clique_search_state& state,
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets,
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

template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
inline bool init(ancestor_clique_search_state& state,
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets,
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

template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
inline bool has_descendant(
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets,
		unsigned int root, unsigned int vertex)
{
	return sets.sets[root].descendants.contains(vertex);
}

template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
bool expand_clique_search_state(
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets,
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

template<bool RecurseChildren, bool TestCompletion, bool ReturnOnCompletion,
	typename StateData, typename BuiltInConstants, typename Formula, typename ProofCalculus, typename... StateArgs>
bool process_clique_search_state(
		search_queue<StateData>& queue,
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets,
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
template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
bool find_largest_disjoint_subset_clique(
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets,
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

template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
inline bool add_non_ancestor_neighbors(
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets,
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
	template<typename BuiltInConstants, typename Formula, typename ProofCalculus, typename Function>
	bool for_each_parent(const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets, unsigned int set, Function apply) const {
		for (const auto& entry : sets.extensional_graph.vertices[set].parents)
			if (!apply(entry.key)) return false;
		for (unsigned int parent : sets.intensional_graph.vertices[set].parents)
			if (!apply(parent)) return false;
		return true;
	}

	template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
	inline unsigned int count(const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets, unsigned int set) const {
		return sets.extensional_graph.vertices[set].parents.size
			 + sets.intensional_graph.vertices[set].parents.length;
	}
};

struct singleton_parent {
	unsigned int parent;

	template<typename BuiltInConstants, typename Formula, typename ProofCalculus, typename Function>
	bool for_each_parent(const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets, unsigned int set, Function apply) const {
		return apply(parent);
	}

	template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
	constexpr unsigned int count(const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets, unsigned int set) const {
		return 1;
	}
};

template<typename BuiltInConstants, typename Formula, typename ProofCalculus, typename SetParents>
bool get_non_ancestor_neighborhoods(
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets,
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

template<typename BuiltInConstants, typename Formula, typename ProofCalculus, typename SetParents>
bool find_largest_disjoint_clique_with_set(
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets,
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
			int priority = sets.sets[next.state->ancestor].set_size;
			for (unsigned int i = 0; i < completed_clique_count; i++)
				priority -= sets.sets[completed_clique[i]].set_size;
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

template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
inline bool find_largest_disjoint_clique_with_set(
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets,
		unsigned int set, unsigned int*& clique, unsigned int& clique_count,
		unsigned int& ancestor_of_clique, int min_priority)
{
	default_parents parents;
	return find_largest_disjoint_clique_with_set(sets, set, parents, clique, clique_count, ancestor_of_clique, min_priority);
}

template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
inline bool find_largest_disjoint_clique_with_set(
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets,
		unsigned int set, unsigned int parent,
		unsigned int*& clique, unsigned int& clique_count,
		unsigned int& ancestor_of_clique, int min_priority)
{
	singleton_parent set_parent = { parent };
	return find_largest_disjoint_clique_with_set(sets, set, set_parent, clique, clique_count, ancestor_of_clique, min_priority);
}

template<typename BuiltInConstants, typename Formula, typename ProofCalculus, typename FirstSetParents, typename SecondSetParents>
bool find_largest_disjoint_clique_with_edge(
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets,
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

template<typename BuiltInConstants, typename Formula, typename ProofCalculus>
inline bool find_largest_disjoint_clique_with_edge(
		const set_reasoning<BuiltInConstants, Formula, ProofCalculus>& sets,
		unsigned int first, unsigned int second,
		unsigned int*& clique, unsigned int& clique_count,
		unsigned int& ancestor_of_clique, int min_priority,
		int stop_priority)
{
	default_parents parents;
	return find_largest_disjoint_clique_with_edge(sets, first, second, parents, parents, clique, clique_count, ancestor_of_clique, min_priority, stop_priority);
}

#endif /* SET_GRAPH_H_ */
