#pragma once

struct MSize {
	unsigned int nrow;
	unsigned int ncol;

	MSize() :nrow(0), ncol(0) {};
	MSize(unsigned int nr, unsigned int nc) : nrow(nr), ncol(nc) {};
	bool IsSquare(){
	    return (nrow==ncol);
	}
	operator bool() {
		return (nrow&&ncol);
	}
	operator int() {
		return nrow * ncol;
	}
	MSize& operator=(const MSize& other) {
		nrow = other.nrow;
		ncol = other.ncol;
		return *this;
	}
};
inline bool operator==(const MSize& l, const MSize& r) { return (l.ncol == r.ncol && l.nrow == r.nrow); }
inline bool operator!=(const MSize& l, const MSize& r) { return (l.ncol != r.ncol || l.nrow != r.nrow); }

inline bool operator<(const MSize& l, const MSize& r) { return (l.ncol*l.nrow < r.ncol*r.nrow); }
inline bool operator>(const MSize& l, const MSize& r) { return (l.ncol*l.nrow > r.ncol*r.nrow); }
inline bool operator<=(const MSize& l, const MSize& r) { return !(l.ncol*l.nrow > r.ncol*r.nrow); }
inline bool operator>=(const MSize& l, const MSize& r) { return !(l.ncol*l.nrow < r.ncol*r.nrow); }