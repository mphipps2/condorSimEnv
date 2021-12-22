/** @file XMLSettingsReader
 *  @brief Function prototypes for XMLSettingsReader
 *
 *  @author Riccardo Longo
 *  @bug No known bugs.
 */

#ifndef XMLSETTINGSREADER_H
#define XMLSETTINGSREADER_H

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

using namespace xercesc;

class XMLSettingsReader {

 public :
      XMLSettingsReader( void );
      virtual ~XMLSettingsReader( void );

      int getBaseNodeCount(std::string _nodeName);
      int getNodeChildCount(std::string _nodeName, int _nodeNumber, std::string _childName);

      std::string getChildValue(std::string _baseNode, int _baseNumber, std::string _childNode);

      void getChildValue(std::string _baseNode, int _baseNumber, std::string _childNode, std::string &_retVal);
      void getChildValue(std::string _baseNode, int _baseNumber, std::string _childNode, bool& _retVal);
      void getChildValue(std::string _baseNode, int _baseNumber, std::string _childNode, double& _retVal);
      void getChildValue(std::string _baseNode, int _baseNumber, std::string _childNode, int& _retVal);

      bool parseFile(std::string _fileName);

 private :
      XercesDOMParser* m_DOMParser;
      DOMElement* m_rootNode;
};

#endif
