#ifndef CONTEXT_H_
#define CONTEXT_H_

#include "lf_utils.h"

struct anaphora_binding {
	unsigned int* indices;

	~anaphora_binding() { free_helper(); }

	static inline void free(anaphora_binding& binding) {
		binding.free_helper();
	}

	static inline bool clone(
			const anaphora_binding& src,
			anaphora_binding& dst,
			unsigned int length)
	{
		if (!dst.init(length)) return false;
		for (unsigned int i = 0; i < length; i++)
			dst.indices[i] = src.indices[i];
		return true;
	}

private:
	inline void free_helper() {
		core::free(indices);
	}

	inline bool init(unsigned int length) {
		indices = (unsigned int*) malloc(max((size_t) 1, sizeof(unsigned int) * length));
		if (indices == nullptr) {
			fprintf(stderr, "anaphora_binding.init ERROR: Out of memory.\n");
			return false;
		}
		return true;
	}

	friend bool init(anaphora_binding&, unsigned int);
};

bool init(anaphora_binding& binding, unsigned int length) {
	return binding.init(length);
}

template<typename Stream>
bool read(anaphora_binding& binding, Stream& in, unsigned int length)
{
	if (!binding.init(length)) {
		return false;
	} else if (!read(binding.indices, in, length)) {
		core::free(binding);
		return false;
	}
	return true;
}

template<typename Stream>
bool write(const anaphora_binding& binding, Stream& out, unsigned int length)
{
	return write(binding.indices, out, length);
}

struct context {
	array<referent_iterator> referent_iterators;
	array<anaphora_binding> bindings;

	context() : referent_iterators(8), bindings(8) { }

	~context() { free_helper(); }

	static inline void free(context& ctx) {
		ctx.free_helper();
		core::free(ctx.referent_iterators);
		core::free(ctx.bindings);
	}

	static inline bool clone(const context& src, context& dst,
			hash_map<const hol_term*, hol_term*>& formula_map)
	{
		if (!array_init(dst.referent_iterators, src.referent_iterators.capacity)) {
			return false;
		} else if (!array_init(dst.bindings, src.bindings.capacity)) {
			core::free(dst.referent_iterators);
			return false;
		}

		for (const referent_iterator& ref_iterator : src.referent_iterators) {
			if (!referent_iterator::clone(ref_iterator, dst.referent_iterators[dst.referent_iterators.length], formula_map)) {
				core::free(dst);
				return false;
			}
			dst.referent_iterators.length++;
		} for (unsigned int i = 0; i < src.bindings.length; i++) {
			const anaphora_binding& binding = src.bindings[i];
			if (!anaphora_binding::clone(binding, dst.bindings[dst.bindings.length], src.referent_iterators[i].anaphora.size)) {
				core::free(dst);
				return false;
			}
			dst.bindings.length++;
		}
		return true;
	}

private:
	inline void free_helper() {
		for (referent_iterator& iter : referent_iterators)
			core::free(iter);
		for (anaphora_binding& binding : bindings)
			core::free(binding);
	}
};

template<typename Formula>
bool get_formula_map(const context& ctx, hash_map<const Formula*, unsigned int>& formula_map)
{
	for (const referent_iterator& iterator : ctx.referent_iterators)
		if (!get_formula_map(iterator, formula_map)) return false;
	return true;
}

template<typename Stream, typename Formula>
bool read(context& ctx, Stream& in, Formula** formulas)
{
	size_t ref_iterator_count;
	if (!read(ref_iterator_count, in)
	 || !array_init(ctx.referent_iterators, ((size_t) 1) << (core::log2(ref_iterator_count == 0 ? 1 : ref_iterator_count) + 1)))
	{
		return false;
	} else if (!array_init(ctx.bindings, ctx.referent_iterators.capacity)) {
		core::free(ctx.referent_iterators);
		return false;
	}

	for (size_t i = 0; i < ref_iterator_count; i++) {
		if (!read(ctx.referent_iterators[ctx.referent_iterators.length], in, formulas)) {
			core::free(ctx);
			return false;
		}
		ctx.referent_iterators.length++;
	} for (size_t i = 0; i < ref_iterator_count; i++) {
		if (!read(ctx.bindings[ctx.bindings.length], in, ctx.referent_iterators[i].anaphora.size)) {
			core::free(ctx);
			return false;
		}
		ctx.bindings.length++;
	}
	return true;
}

template<typename Stream, typename Formula>
bool write(const context& ctx, Stream& out,
		const hash_map<const Formula*, unsigned int>& formula_map)
{
	if (!write(ctx.referent_iterators.length, out))
		return false;

	for (size_t i = 0; i < ctx.referent_iterators.length; i++) {
		if (!write(ctx.referent_iterators[i], out, formula_map))
			return false;
	} for (size_t i = 0; i < ctx.bindings.length; i++) {
		if (!write(ctx.bindings[i], out, ctx.referent_iterators[i].anaphora.size))
			return false;
	}
	return true;
}

#endif /* CONTEXT_H_ */
