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

	static inline bool clone(
			const context& src, context& dst)
	{
		if (!array_init(dst.referent_iterators, src.referent_iterators.capacity)) {
			return false;
		} else if (!array_init(dst.bindings, src.bindings.capacity)) {
			core::free(dst.referent_iterators);
			return false;
		}

		for (const referent_iterator& ref_iterator : src.referent_iterators) {
			if (!referent_iterator::clone(ref_iterator, dst.referent_iterators[dst.referent_iterators.length])) {
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

#endif /* CONTEXT_H_ */
