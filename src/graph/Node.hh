#ifndef NODE_H
#define NODE_H

#include <vector>
#include <string>
#include <memory>

using namespace std;

namespace bnc {

    /**
     * Define the property on Edges
     * e.g. for edge "X~norm(Y,sd)" ("Y->X"), the property
     * should be 1 for X and Y because Y is the first parameter
     * (mu) of the distribution of X;
     */
    typedef int ParamIndex;

    template<class DistType=string, class EdgeProp=ParamIndex>
    class Node {
    public:
	Node(const DistType ) {};
	~Node() {};

	const vector<share_ptr<Node>>& inNodes() const {
	    return _inNodes;
	}

	const vector<share_ptr<Node>>& outNodes() const {
	    return _outNodes;
	}

	const vector<share_ptr<DistType>>& inTypes() const {
	    return _inTypes;
	}

	const vector<share_ptr<DistType>>& outTypes() const {
	    return _outTypes
	}

	const vector<int>& inParamIndex() const {
	    return _inParamIndex;
	}

	const vector<int>& outParamIndex() const {
	    return _outParamIndex;
	}

    private:
	vector<share_ptr<Node>> _inNodes;
	vector<share_ptr<Node>> _outNodes;
	vector<share_ptr<DistType>> _inTypes;
	vector<share_ptr<DistType>> _outTypes;
	vector<int> _inParamIndex;
	vector<int> _outParamIndex;

    };

    inline bool isRoot(const share_ptr<Node>& node) {
	return node->inNodes().size() == 0;
    }

    inline bool isLeaf(const share_ptr<Node>& node) {
	return node->outNodes().size() == 0;
    }
    
} // namespace bnc

#endif /* NODE_H */
