// (c) Copyright Juergen Hunold 2016
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MODULE QtDataVisualization

#include <QtDataVisualization>

#include <QGuiApplication>

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE (defines)
{
    BOOST_CHECK_EQUAL(BOOST_IS_DEFINED(QT_CORE_LIB), true);
    BOOST_CHECK_EQUAL(BOOST_IS_DEFINED(QT_GUI_LIB), true);
    BOOST_CHECK_EQUAL(BOOST_IS_DEFINED(QT_DATAVISUALIZATION_LIB), true);
}

BOOST_AUTO_TEST_CASE( datavisualization )
{
    QGuiApplication app(boost::unit_test::framework::master_test_suite().argc,
                        boost::unit_test::framework::master_test_suite().argv);

    QtDataVisualization::Q3DBars graph;

    graph.setShadowQuality(QtDataVisualization::QAbstract3DGraph::ShadowQualitySoftMedium);
    graph.activeTheme()->setBackgroundEnabled(false);
    graph.activeTheme()->setLabelBackgroundEnabled(true);
}
