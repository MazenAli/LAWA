#include <iostream>
#include <htucker/htucker.h>

typedef double                          T;
typedef htucker::HTuckerTree<T>         HTTree;
typedef htucker::HTuckerTreeNode<T>     HNode;
typedef htucker::GeneralTreeIterator
        <HNode>                         TreeIt;
typedef htucker::GeneralTreeNode<
        htucker::HTuckerTreeNode<T> >   GNode;

int
main()
{
    int     d = 8;
    HTTree  X(d);
    GNode   *node;

    std::cout << "--- Here begins the tree iteration ---\n";
    std::cout << "The tree is =>\n";
    X.print();
    std::cout << std::endl;
    for(TreeIt tit=X.getGeneralTree().end();
               tit>=X.getGeneralTree().begin(); tit--) {
        node  = tit.getNode();
        std::cout << "Current node =>\n";
        node->printnode();
        std::cout << "\nfirst child =>\n";
        if (node->getfirstChild())
        node->getfirstChild()->printnode();
        std::cout << "\nlast child =>\n";
        if (node->getlastChild())
        node->getlastChild()->printnode();
        std::cout << "\nprevious sibling=>\n";
        if (node->getpreviousSibling())
        node->getpreviousSibling()->printnode();
        std::cout << "\nnext sibling=>\n";
        if (node->getnextSibling())
        node->getnextSibling()->printnode();
        std::cout << "\nlevel left=>\n";
        if (node->getlevelleft())
        node->getlevelleft()->printnode();
        std::cout << "\nlevel right=>\n";
        if (node->getlevelright())
        node->getlevelright()->printnode();
        std::cout << "\n|\n|\n|\n|\n";
    }

    return 0;
}
