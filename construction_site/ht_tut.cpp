#include <iostream>
#include <htucker/htucker.h>
#include <flens/flens.cxx>

typedef double                                                  T;
typedef flens::DenseVector<flens::Array<T> >                    DenseVector;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
                                                                Matrix;


int
main()
{
    htucker::HTuckerTree<T>         X(5);
    htucker::HTuckerTree<T>         _X;
    htucker::HTuckerTree<T>         diff;
    DenseVector                     b(3);
    DenseVector                     _i(5);
    htucker::DimensionIndex         i(5);

    b = 1, 2, 3;
    htucker::DenseVectorList<T>     list(b);
    b = 1, 2, 3.0000001;
    list.add(b);
    b = 3, 4, 5;
    list.add(b);
    b = 8, 9, 10;
    list.add(b);
    b = 3, 4, 3;
    list.add(b);
    b = 6, 8, 5.9999999;
    list.add(b);
    b = 8, 9, 10;
    list.add(b);
    b = 3, 4, 5;
    list.add(b);

    b = 2, 1, 5;
    list.add(b);
    b = 4, 2, 9.999995;
    list.add(b);


    X.generateTofElementary(list, 2, 5);
    _i = 1, 1, 1, 1, 1;
    i.setValue(_i);
    std::cout << i << " =>\n";
    std::cout << X.evaluate(i) << std::endl;
    _i = 3, 3, 3, 3, 3;
    i.setValue(_i);
    std::cout << i << " =>\n";
    std::cout << X.evaluate(i) << std::endl;

    int min = 1;
    int max = 3;
    i.setRandom(min, max);
    std::cout << i << " =>\n";
    std::cout << X.evaluate(i) << std::endl;

    i.setRandom(min, max);
    std::cout << i << " =>\n";
    std::cout << X.evaluate(i) << std::endl;

    i.setRandom(min, max);
    std::cout << i << " =>\n";
    std::cout << X.evaluate(i) << std::endl;

    std::cout << "\n printing frames\n";
    X.print_w_UorB();
    X.print_info();
    std::cout << "\n printing values\n";
    X.print_values();

    htucker::GeneralTreeNode<htucker::HTuckerTreeNode<T> >* node;
    htucker::GeneralTreeNode<htucker::HTuckerTreeNode<T> >* save;
    for (htucker::GeneralTreeIterator<htucker::HTuckerTreeNode<T> >
         tit = X.getGeneralTree().end();
         tit >= X.getGeneralTree().begin(); tit--) {
        save = tit.getNode();
        std::cout << "This is node =>\n";
        save->printnode();
        std::cout << "With level =>\n";
        std::cout << save->level() << std::endl;

        if (save->isLeaf()) {
            Matrix      M = save->getContent()->getUorB();
            Matrix      U, V;
            DenseVector s;
            flens::svd(M, s, U, V);
            std::cout << "The sv's are =>\n";
            std::cout << s << std::endl;
            std::cout << "Corresponding basis U =>\n";
            std::cout << U << std::endl;
        }


        std::cout << "\nIt's parents is =>\n";
        node = save->getParent();
        if (node == NULL) {
            std::cout << "NULL\n";
        } else {
            node->printnode();
        }
        std::cout << "\nIt's first child is =>\n";
        node = save->getfirstChild();
        if (node == NULL) {
            std::cout << "NULL\n";
        } else {
            node->printnode();
        }

        std::cout << "\nIt's last child is =>\n";
        node = save->getlastChild();
        if (node == NULL) {
            std::cout << "NULL\n";
        } else {
            node->printnode();
        }

        std::cout << "\nIt's last sibling is =>\n";
        node = save->getpreviousSibling();
        if (node == NULL) {
            std::cout << "NULL\n";
        } else {
            node->printnode();
        }

        std::cout << "\nIt's next sibling is =>\n";
        node = save->getnextSibling();
        if (node == NULL) {
            std::cout << "NULL\n";
        } else {
            node->printnode();
        }

        std::cout << "\nIt's left level is =>\n";
        node = save->getlevelleft();
        if (node == NULL) {
            std::cout << "NULL\n";
        } else {
            node->printnode();
        }

        std::cout << "\nIt's right level is =>\n";
        node = save->getlevelright();
        if (node == NULL) {
            std::cout << "NULL\n";
        } else {
            node->printnode();
        }
    }


    X.orthogonalize();
    std::cout << "\n printing values\n";
    X.print_values();
    std::cout << "\nmax rank\n";
    std::cout << X.max_rank() << std::endl;
    std::cout << "\naverage rank\n";
    std::cout << X.average_rank() << std::endl;
    std::cout << "\neffective rank\n";
    std::cout << X.effective_rank() << std::endl;

    std::cout << "Post truncation with orthogonalization to rank 1 =>\n";
//    _X.print_w_UorB();
//    X.print_values();
//    _X.print_values();

    std::cout << "In place truncation =>\n";
//    X.print_w_UorB();
//    _X = htucker::truncate(X, 1);
//    X.truncate(1, true);
    _X = X;
    X.print_values();
    double eps = 1e-03;
    _X.truncate(eps, true);
    _X.print_values();
//    X.print_w_UorB();
    _X.print_w_UorB();

    diff = X - _X;
    diff.orthogonalize();
    std::cout << "L2 error = " << diff.L2normorthogonal() << std::endl;

    /*
    std::cout << "\n---Sorting---\n";
    DenseVector s(7);
    s = 1, 4, 5, 10, 8, 5, 3;
    std::cout << "Old s =>\n" << s << std::endl;
    flens::DenseVector<flens::Array<int> > rho;
    flens::DenseVector<flens::Array<int> > blub(5);
    flens::sort(s, rho);
    std::cout << "New s =>\n" << s << std::endl;
    std::cout << "rho =>\n" << rho << std::endl;

    ++blub(3);
    std::cout << "blub=>\n" << blub << std::endl;
    */

    return 0;
}
