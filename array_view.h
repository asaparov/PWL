#ifndef ARRAY_VIEW_H_
#define ARRAY_VIEW_H_

namespace detail {
	template<typename A, typename C> static auto test_index_operator(int32_t) ->
			decltype(C(std::declval<A>()[0]), std::true_type{});
	template<typename A, typename C> static auto test_index_operator(int64_t) -> std::false_type;
}

template<typename T, typename ReturnType> struct has_index_operator : decltype(::detail::test_index_operator<T, ReturnType>(0))::type {};

template<bool Unique, typename T>
inline void add_sorted(array<T>& list, const T& element) {
	unsigned int index = linear_search(list.data, element, 0, list.length);
	if (Unique && index < list.length && list[index] == element) return;
	shift_right(list.data, list.length, index);
	list[index] = element;
	list.length++;
}

template<typename T>
struct array_view {
	T* array;
	unsigned int length;

	array_view(T* array, unsigned int length) : array(array), length(length) { }

	inline T& operator[] (size_t index) {
		return array[index];
	}

	inline const T& operator[] (size_t index) const {
		return array[index];
	}

	inline const T* begin() const {
		return array;
	}

	inline const T* end() const {
		return array + length;
	}

	inline unsigned int size() const {
		return length;
	}
};

template<typename T>
array_view<T> make_array_view(T* array, unsigned int length) {
	return array_view<T>(array, length);
}

template<typename T>
struct indexed_array_view {
	T* array;
	const unsigned int* indices;
	unsigned int length;

	indexed_array_view(T* array, const unsigned int* indices, unsigned int length) : array(array), indices(indices), length(length) { }

	inline T& operator[] (size_t index) {
		return array[indices[index]];
	}

	inline const T& operator[] (size_t index) const {
		return array[indices[index]];
	}

	inline unsigned int size() const {
		return length;
	}
};

template<typename T>
indexed_array_view<T> make_indexed_array_view(T* array, unsigned int* indices, unsigned int length) {
	return indexed_array_view<T>(array, indices, length);
}

template<typename Array, typename T>
struct lookup_table_array_view {
	static_assert(has_index_operator<Array, T>::value, "`Array` does not have an index operator that returns type `T`");

	Array* arrays;
	const unsigned int* indices;
	unsigned int length;

	lookup_table_array_view(Array* arrays, const unsigned int* indices, unsigned int length) : arrays(arrays), indices(indices), length(length) { }

	inline T& operator[] (size_t index) {
		return arrays[index][indices[index]];
	}

	inline const T& operator[] (size_t index) const {
		return arrays[index][indices[index]];
	}

	inline unsigned int size() const {
		return length;
	}
};

template<typename T, template<typename> class Array>
struct prepended_array_view {
	T& first;
	const Array<T>& second;

	prepended_array_view(T& first, const Array<T>& second) : first(first), second(second) { }

	inline T& operator[] (size_t index) {
		if (index == 0) return first;
		else return second[index - 1];
	}

	inline const T& operator[] (size_t index) const {
		if (index == 0) return first;
		else return second[index - 1];
	}

	inline unsigned int size() const {
		return 1 + second.size();
	}
};

template<typename T, template<typename> class Array>
prepended_array_view<T, Array> make_prepended_array_view(T& first, const Array<T>& second) {
	return prepended_array_view<T, Array>(first, second);
}

template<typename T, template<typename> class Array>
struct appended_array_view {
	const Array<T>& first;
	T& second;

	appended_array_view(const Array<T>& first, T& second) : first(first), second(second) { }

	inline T& operator[] (size_t index) {
		if (index == first.size()) return second;
		else return first[index];
	}

	inline const T& operator[] (size_t index) const {
		if (index == first.size()) return second;
		else return first[index];
	}

	inline unsigned int size() const {
		return first.size() + 1;
	}
};

template<typename T, template<typename> class Array>
appended_array_view<T, Array> make_appended_array_view(const Array<T>& first, T& second) {
	return appended_array_view<T, Array>(first, second);
}

template<typename T>
struct excluded_array_view {
	T* elements;
	unsigned int length;
	unsigned int excluded_index;

	excluded_array_view(T* elements, unsigned int original_length, unsigned int excluded_index) :
			elements(elements), length(original_length - 1), excluded_index(excluded_index) { }

	inline T& operator[] (size_t index) {
		if (index < excluded_index)
			return elements[index];
		else return elements[index + 1];
	}

	inline const T& operator[] (size_t index) const {
		if (index < excluded_index)
			return elements[index];
		else return elements[index + 1];
	}

	inline unsigned int size() const {
		return length;
	}
};

template<typename T>
inline excluded_array_view<T> make_excluded_array_view(T* elements, unsigned int original_length, unsigned int excluded_index) {
	return excluded_array_view<T>(elements, original_length, excluded_index);
}

template<typename T>
struct repeated_array_view {
	T& repeated_element;
	unsigned int length;

	repeated_array_view(T& repeated_element, unsigned int length) :
			repeated_element(repeated_element), length(length) { }

	inline T& operator[] (size_t index) {
		return repeated_element;
	}

	inline const T& operator[] (size_t index) const {
		return repeated_element;
	}

	inline unsigned int size() const {
		return length;
	}
};

template<typename T>
inline repeated_array_view<T> make_repeated_array_view(T& repeated_element, unsigned int length) {
	return repeated_array_view<T>(repeated_element, length);
}

#endif /* ARRAY_VIEW_H_ */
