#include "control/diagnostic_tool/solver_structure_visualizer.h"

#include "data_management/output_connector_data.h"
#include "output_connector_data_transfer/output_connection.h"
#include "output_writer/generic.h"
#include "utility/string_utility.h"

const int VARIABLES_LINE_LENGTH = 45;  // number of characters for the solver structure and variables, afterwards there will be connection lines 

//! iterate over all nested solvers
void SolverStructureVisualizer::DiagramGenerator::
generateDiagramRecursion(std::stringstream &result, std::vector<std::vector<int>> &internalConnectionLinesCurrentSolver, std::vector<int> &slotLineNos,
                         int depth, bool isFirstChildSolver, bool isLastChildSolver)
{
  if (!currentSolver_)
    return;

  // prepare the output connections
  SolverStructureVisualizer::parseOutputConnection(currentSolver_);

  slotLineNos.clear();

  std::stringstream lineStart;
  for (int i = 0; i < depth; i++)
    lineStart << "│ ";

  result << lineStart.str();
  if (isLastChildSolver && currentSolver_->outputSlots.empty() && currentSolver_->outputConnections.empty())
    result << "└";
  else
    result << "├";
    //result << "┌";
  result << "── " << currentSolver_->name << "\n";

  if (currentSolver_->description != "")
  {
    result << lineStart.str() << "│   (\"" << currentSolver_->description << "\")" << "\n";
  }

  // print the output slots of the solver, store the line nos of each slot to slotLineNos
  if (!currentSolver_->outputSlots.empty())
  {
    result << lineStart.str() << "│  output slots: \n";

    // loop over the output connector slots of the current solver
    for (int i = 0; i < currentSolver_->outputSlots.size(); i++)
    {
      // print field variable and component name of the output slot
      std::stringstream s;
      std::string fieldVariableName = currentSolver_->outputSlots[i].fieldVariableName;
      std::string componentName = currentSolver_->outputSlots[i].componentName;

      // only add component name if it is different from the variable name
      if (fieldVariableName != componentName)
      {
        s << lineStart.str() << "│  " << fieldVariableName;

        // only add component name if there are multiple components or, or the component is not not named "0" which is the standard
        if (currentSolver_->outputSlots[i].nComponents > 1 || componentName != "0")
          s << "." << componentName;
      }
      else
      {
        s << lineStart.str() << "│  " << fieldVariableName;
      }

      // add note about variableNo
      if (currentSolver_->outputSlots[i].variableNo != 1)
        s << " (in variable" << currentSolver_->outputSlots[i].variableNo << ")";

      std::string outputSlotString = s.str();

      int currentLineLength = StringUtility::stringLength(outputSlotString)+2;
      const int requiredLineLength = VARIABLES_LINE_LENGTH;

      // shorten string if neccessary
      if (currentLineLength >= requiredLineLength)
      {
        outputSlotString = outputSlotString.substr(0, outputSlotString.length() - (currentLineLength - requiredLineLength));
        currentLineLength = StringUtility::stringLength(outputSlotString)+2;
        s.str("");
        s << outputSlotString;
      }

      s << std::string(requiredLineLength - currentLineLength, ' ')
        << "── ¤" << i << "\n";
      result << s.str();

      // store line no
      std::string currentResult = result.str();
      int currentLineNo = std::count(currentResult.begin(), currentResult.end(), '\n')-1;

      //LOG(DEBUG) << "current Result [" << currentResult << "] contains " << currentLineNo << " line breaks";
      slotLineNos.push_back(currentLineNo);
    }

    result << lineStart.str() << "│\n";
  }

  // maybe todo for later: show also internal connections of all slots to the first subsolver

  // print the output slot connections of the solver
  std::vector<solver_t::OutputConnectionRepresentation> &outputConnections = currentSolver_->outputConnections;
  if (!outputConnections.empty())
  {
    result << lineStart.str() << "│  slot connections: \n";
    for (int i = 0; i < outputConnections.size(); i++)
    {
      if (outputConnections[i].type == solver_t::OutputConnectionRepresentation::output_connection_t::ba)
      {
        result << lineStart.str() << "│  " <<  outputConnections[i].toSlot << "¤ <─ ¤" << outputConnections[i].fromSlot << "\n";
      }
      else
      {
        result << lineStart.str() << "│  " <<  outputConnections[i].fromSlot << "¤";

        switch (outputConnections[i].type)
        {
        case solver_t::OutputConnectionRepresentation::output_connection_t::ba:
        case solver_t::OutputConnectionRepresentation::output_connection_t::ab:
          result << " ─> ";
          break;
        case solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalReuse:
          result << " <═> ";  // or use ⇔ ?
          break;
        case solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalCopy:
          result << " <─> ";
          break;

        };

        result << "¤" << outputConnections[i].toSlot << "\n";
      }
    }
    result << lineStart.str() << "│\n";
  }

  // get the line nos of the slots for nested solvers 0 and 1
  std::vector<std::shared_ptr<solver_t>> &children = currentSolver_->children;
  std::vector<std::vector<int>> slotLineNosTerm(children.size());
  std::vector<std::vector<std::vector<int>>> internalConnectionLinesCurrentSolverTerm(children.size());

  std::shared_ptr<solver_t> ownSolver = currentSolver_;

  int i = 0;
  // loop over child solvers
  for (std::vector<std::shared_ptr<solver_t>>::iterator iter = children.begin();
       iter != children.end(); iter++, i++)
  {
    currentSolver_ = (*iter);

    bool isFirstChildSolver = i == 0;
    bool isLastChildSolver = i == children.size()-1;

    // recursively call generateDiagram on subsolvers
    generateDiagramRecursion(result, internalConnectionLinesCurrentSolverTerm[i], slotLineNosTerm[i], depth+1, isFirstChildSolver, isLastChildSolver);

    VLOG(1) << std::string(depth*2, ' ') << "child " << i << ", i:" << internalConnectionLinesCurrentSolverTerm[i];

    // add end of subsolver
    if (isLastChildSolver)
      result << lineStart.str() << "└\n";
    else
      result << lineStart.str() << "│\n";

    // internal connection lines of the the 2nd and further solvers are finished, the one of the first are passed on to the parent solver
    if (isFirstChildSolver && ownSolver->hasInternalConnectionToFirstNestedSolver)
    {
      VLOG(1) << std::string(depth*2, ' ') << "(" << ownSolver->name << ") first child solver gives " << internalConnectionLinesCurrentSolverTerm[i] << " own are: " << internalConnectionLinesCurrentSolver;
      internalConnectionLinesCurrentSolver.insert(internalConnectionLinesCurrentSolver.end(), internalConnectionLinesCurrentSolverTerm[i].begin(), internalConnectionLinesCurrentSolverTerm[i].end());
      VLOG(1) << std::string(depth*2, ' ') << "now own: " << internalConnectionLinesCurrentSolver;
    }
    else 
    {
      VLOG(1) << std::string(depth*2, ' ') << "(" << ownSolver->name << ") " << i << ". child solver gives " << internalConnectionLinesCurrentSolverTerm[i] << " all are: " << internalConnectionLinesAll_;
      internalConnectionLinesAll_.insert(internalConnectionLinesAll_.end(), internalConnectionLinesCurrentSolverTerm[i].begin(), internalConnectionLinesCurrentSolverTerm[i].end());
      VLOG(1) << std::string(depth*2, ' ') << "now all: " << internalConnectionLinesAll_;
    }
  }

  // loop over the output connector slots of the current solver
  for (int i = 0; i < ownSolver->outputSlots.size(); i++)
  {
    // add slot to internal connection lines
    if (i >= internalConnectionLinesCurrentSolver.size())
    {
      internalConnectionLinesCurrentSolver.resize(i+1);
    }
    internalConnectionLinesCurrentSolver[i].push_back(slotLineNos[i]);
    VLOG(1) << std::string(depth*2, ' ') << "(" << ownSolver->name << ") add slot " << slotLineNos[i] << ", now i: " << internalConnectionLinesCurrentSolver;
  }

  VLOG(1) << std::string(depth*2, ' ') << "internalConnectionLinesCurrentSolver:" << internalConnectionLinesCurrentSolver;
  VLOG(1) << std::string(depth*2, ' ') << "internalConnectionLinesAll_:" << internalConnectionLinesAll_;

  // add connections between slots of nested solvers to externalConnectionLines
  for (int i = 0; i < outputConnections.size(); i++)
  {
    if (slotLineNosTerm.size() >= 2)
    {
      if ((outputConnections[i].type == solver_t::OutputConnectionRepresentation::output_connection_t::ab
        || outputConnections[i].type == solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalCopy
        || outputConnections[i].type == solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalReuse)
        && slotLineNosTerm[0].size() > outputConnections[i].fromSlot
        && slotLineNosTerm[1].size() > outputConnections[i].toSlot
      )
      {
        int lineNoFrom = slotLineNosTerm[0][outputConnections[i].fromSlot];
        int lineNoTo = slotLineNosTerm[1][outputConnections[i].toSlot];
        solver_t::OutputConnectionRepresentation::output_connection_t connectionType = outputConnections[i].type;

        externalConnectionLines_.push_back(std::make_tuple(
          lineNoFrom, lineNoTo, 0, connectionType
        ));
      }
      else if (outputConnections[i].type == solver_t::OutputConnectionRepresentation::output_connection_t::ba
        && slotLineNosTerm[1].size() > outputConnections[i].fromSlot
        && slotLineNosTerm[0].size() > outputConnections[i].toSlot
      )
      {
        int lineNoFrom = slotLineNosTerm[0][outputConnections[i].toSlot];
        int lineNoTo = slotLineNosTerm[1][outputConnections[i].fromSlot];
        solver_t::OutputConnectionRepresentation::output_connection_t connectionType = outputConnections[i].type;

        externalConnectionLines_.push_back(std::make_tuple(
          lineNoFrom, lineNoTo, 0, connectionType
        ));
      }
      else
      {
        VLOG(1) << "slotLineNosTerm is not good. fromSlot: " << outputConnections[i].fromSlot
          << ", toSlot: " << outputConnections[i].toSlot << ", slot line nos term0: " << slotLineNosTerm[0]
          << ", term1: " << slotLineNosTerm[1];
      }
    }
    else
    {
      VLOG(1) << "slotLineNosTerm is too short: " << slotLineNosTerm.size() << ", but there are " << outputConnections.size() << " outputConnections.";
    }
  }

  // if the own solver does not have output slots but it has only a single child, use their output slots,
  // this is true for MultipleInstances
  if (slotLineNos.empty() && children.size() == 1)
  {
    slotLineNos = slotLineNosTerm[0];
  }
  slotLineNosAll_.insert(slotLineNosAll_.end(), slotLineNos.begin(), slotLineNos.end());
}

//! initialize the diagram generator, afterwards, generateDiagram can be called
void SolverStructureVisualizer::DiagramGenerator::
initialize(std::shared_ptr<solver_t> solverRoot)
{
  currentSolver_ = solverRoot;
}

std::string SolverStructureVisualizer::DiagramGenerator::
generateDiagramWithoutConnectionLines()
{
  std::stringstream diagramWithoutConnectionLines;

  std::vector<std::vector<int>> internalConnectionLinesCurrentSolver;
  std::vector<int> slotLineNos;

  // call the implementation 
  generateDiagramRecursion(diagramWithoutConnectionLines, internalConnectionLinesCurrentSolver, slotLineNos, 0, false, false);

  // add the rest of the internal connection lines that are still storade under internalConnectionLinesCurrentSolverTerm to internalConnectionLinesAll
  internalConnectionLinesAll_.insert(internalConnectionLinesAll_.end(), internalConnectionLinesCurrentSolver.begin(), internalConnectionLinesCurrentSolver.end());

  return std::string(diagramWithoutConnectionLines.str());
}

std::string SolverStructureVisualizer::DiagramGenerator::
generateFinalDiagram()
{
  std::string diagram;

  //! externalConnectionLines_ contains the following information: <lineNoFrom, lineNoTo, lineColumn, lineType>
  // lineColumn is the horizontal position of the vertical data connection line

  // call generateDiagram on the nested solvers, this also fills the variables
  // externalConnectionLines_, internalConnectionLinesAll_ and slotLineNosAll_
  std::string diagramWithoutConnectionLinesString = generateDiagramWithoutConnectionLines();

  // set the "lineColumn" entry in externalConnectionLines
  // loop over externalConnectionLines
  for (int i = 0; i < externalConnectionLines_.size(); i++)
  {
    int startLineNo = std::get<0>(externalConnectionLines_[i]);
    int endLineNo = std::get<1>(externalConnectionLines_[i]);
    int lineColumn = 0;

    // count number of starting lines inside this line
    for (int j = 0; j < externalConnectionLines_.size(); j++)
    {
      if (i == j)
        continue;
      int startLineNo2 = std::get<0>(externalConnectionLines_[j]);
      int endLineNo2 = std::get<1>(externalConnectionLines_[j]);

      if (startLineNo < startLineNo2 && startLineNo2 < endLineNo)
        lineColumn++;
      else if (startLineNo == startLineNo2 && endLineNo2 < endLineNo)
        lineColumn++;
      else if (startLineNo < startLineNo2 && endLineNo2 == endLineNo)
        lineColumn++;
    }
    std::get<2>(externalConnectionLines_[i]) = lineColumn;
  }

  // erase and remove internal connections that consist of one slot only
  internalConnectionLinesAll_.erase(
    std::remove_if(internalConnectionLinesAll_.begin(), internalConnectionLinesAll_.end(), 
    [](std::vector<int> connection){
      return connection.size() == 1;
    }), 
    internalConnectionLinesAll_.end()
  );

  LOG(DEBUG) << "externalConnectionLines_: " << externalConnectionLines_;
  LOG(DEBUG) << "internalConnectionLines: " << internalConnectionLinesAll_;

  // determine number of columns that are needed for internal connections
  int nColumnsForInternalConnectionLines = 0;
  int nLinesDiagram = std::count(diagramWithoutConnectionLinesString.begin(), diagramWithoutConnectionLinesString.end(), '\n');
  LOG(DEBUG) << "nLinesDiagram: " << nLinesDiagram;
    
  std::map<int,int> internalConnectionLineColumns;   //< key = index into internalConnectionLinesAll_, value = columnNo of that line
  for (int currentLineNo = 0; currentLineNo < nLinesDiagram; currentLineNo++)
  {
    // determine active connection lines 
    std::set<int> activeLinesIndices;
    int activeLineIndex = 0;

    // loop over all connection lines that are present in the current row
    for (const std::vector<int> &connectionLine : internalConnectionLinesAll_)
    {
      int lineNoFirst = *std::min_element(connectionLine.begin(), connectionLine.end());
      int lineNoLast = *std::max_element(connectionLine.begin(), connectionLine.end());
      
      // if the currently considered connectionLine is active, i.e. present in the current line "currentLineNo"
      if (lineNoFirst <= currentLineNo && currentLineNo <= lineNoLast)
      {
        // ensure that the connectionLine gets a column assigned if it does not have one already
        
        // if the current activeLine has no assigned column yet
        if (internalConnectionLineColumns.find(activeLineIndex) == internalConnectionLineColumns.end())
        {
          VLOG(1) << "find new column for internal connection line with index " << activeLineIndex << ": " << internalConnectionLinesAll_[activeLineIndex];

          // loop over possible columns until a free column is found that can be assign to the current activeLine
          for (int columnNo = 0; ; columnNo++)
          {
            // check if current column is still free
            bool columnIsFree = true;
            for (std::pair<int,int> columnNoConnectionLine : internalConnectionLineColumns)
            {
              if (columnNoConnectionLine.second == columnNo)
              {
                VLOG(1) << "column " << columnNo << " is occupied by line " << columnNoConnectionLine.first << ": " << internalConnectionLinesAll_[columnNoConnectionLine.first];
                columnIsFree = false;
                break;
              }
            }
            if (columnIsFree)
            {
              VLOG(1) << "found free column " << columnNo;

              // assign column to this connection line
              internalConnectionLineColumns[activeLineIndex] = columnNo;
              nColumnsForInternalConnectionLines = std::max(nColumnsForInternalConnectionLines, columnNo+1);
              break;
            }
          }
        }
      }

      activeLineIndex++;
    }
  }

  LOG(DEBUG) << "internalConnectionLineColumns: " << internalConnectionLineColumns;
  LOG(DEBUG) << "nColumnsForInternalConnectionLines: " << nColumnsForInternalConnectionLines;

  std::size_t pos = 0;
  diagram += "Solver structure: \n\n";

  // print actual diagram, into std::string diagram
  // loop over lines
  for (int currentLineNo = 0; pos < diagramWithoutConnectionLinesString.length(); currentLineNo++)
  {
    // determine current line by iterating over substring between "\n"s
    int nextPos = diagramWithoutConnectionLinesString.find("\n",pos+1);
    std::string line;

    if (pos == 0)
      pos = -1;

    if (nextPos != std::string::npos)
    {
      line = diagramWithoutConnectionLinesString.substr(pos+1, nextPos-pos-1);
    }
    else
    {
      line = diagramWithoutConnectionLinesString.substr(pos+1);
    }

    pos = nextPos;

    // here, line contains the current line that should be handled without "\n"

    // fill current line with spaces such that it has the correct length
    int lineLength = StringUtility::stringLength(line);

    // if the line is too short, add space until is has the correct length
    for (int i = 0; i < VARIABLES_LINE_LENGTH - lineLength; i++)
    {
      line += " ";
    }

#if 1
    // add internal connection lines

    std::string internalConnectionLinesColumns;

    // they are saved in internalConnectionLinesAll_
    // determine number of vertical lines in current row
    for (int currentColumnNo = 0; currentColumnNo < nColumnsForInternalConnectionLines; currentColumnNo++)
    {
      bool connectionLineFoundAtCurrentColumn = false;
      for (std::pair<int,int> internalConnectionLineColumn : internalConnectionLineColumns)
      {
        int columnNo = internalConnectionLineColumn.second;

        if (currentColumnNo == columnNo)
        {
          VLOG(1) << "currentColumnNo " << currentColumnNo << ", connection line " << internalConnectionLineColumn.first << "/" << internalConnectionLinesAll_.size();

          const std::vector<int> &connectionLine = internalConnectionLinesAll_[internalConnectionLineColumn.first];
          
          int lineNoFirst = *std::min_element(connectionLine.begin(), connectionLine.end());
          int lineNoLast = *std::max_element(connectionLine.begin(), connectionLine.end());
          
          VLOG(1) << "line " << currentLineNo << ", connection [" << lineNoFirst << "," << lineNoLast << "]";

          if (lineNoFirst <= currentLineNo && currentLineNo <= lineNoLast)
          {
            // check if it is connected to the current row
            if (std::find(connectionLine.begin(), connectionLine.end(), currentLineNo) != connectionLine.end())
            {
              internalConnectionLinesColumns += "+";
            }
            else
            {
              internalConnectionLinesColumns += ":";
            }

            connectionLineFoundAtCurrentColumn = true;
          }
          break;
        }
      }

      if (!connectionLineFoundAtCurrentColumn)
        internalConnectionLinesColumns += " ";
    }

    std::string firstPart;
    int posFirstPart = 0;
    for(int i = 0; i < line.size(); i++)
    {
      firstPart += line[i];
      if (StringUtility::stringLength(firstPart) == VARIABLES_LINE_LENGTH-2)
      {
        posFirstPart = i;
        break;
      }
    }
    line = firstPart + std::string(" ") + internalConnectionLinesColumns + line.substr(posFirstPart+1);

    //line += internalConnectionLinesColumns;
#endif

    // prepare for external connection lines
    // find out vertial data slot connection lines that go through this output line
    std::vector<std::tuple<int,bool,bool,solver_t::OutputConnectionRepresentation::output_connection_t>> verticalLinesInCurrentRow;    //! <lineColumn,thisLineStartsHere,thisLineEndsHere,type>, lineColumn is the column in rhe current row where the vertical lines passes through
    bool lineStartsHere = false;
    bool lineEndsHere = false;
    bool unconnectedSlotHere = false;
    solver_t::OutputConnectionRepresentation::output_connection_t lineType = solver_t::OutputConnectionRepresentation::output_connection_t::ab;
    int lineColumn = 0;
#if 0
    // get maximum lineColumn
    int maximumLinePosition = 0;
    for (int i = 0; i < externalConnectionLines_.size(); i++)
    {
      maximumLinePosition = std::max((int)(std::get<2>(externalConnectionLines_[i])), maximumLinePosition);
    }
#endif
    // loop over all externalConnectionLines_ and see if they pass through the current line
    // collect all those rows into data structure verticalLinesInCurrentRow
    for (int i = 0; i < externalConnectionLines_.size(); i++)
    {
      // start and end row of current connection line
      int startLineNo = std::get<0>(externalConnectionLines_[i]);
      int endLineNo = std::get<1>(externalConnectionLines_[i]);

      // if the connection line starts or ends in the current row, there will be a horizontal line
      if (currentLineNo == startLineNo)
      {
        lineType = std::get<3>(externalConnectionLines_[i]);
        lineStartsHere = true;
        lineColumn = std::get<2>(externalConnectionLines_[i]);
        verticalLinesInCurrentRow.push_back(std::make_tuple(lineColumn, true, false, lineType));
      }
      else if (currentLineNo == endLineNo)
      {
        lineType = std::get<3>(externalConnectionLines_[i]);
        lineEndsHere = true;
        lineColumn = std::get<2>(externalConnectionLines_[i]);
        verticalLinesInCurrentRow.push_back(std::make_tuple(lineColumn, false, true, lineType));
      }
      else if (startLineNo < currentLineNo && currentLineNo < endLineNo)
      {
        // vertical line passes through current row
        solver_t::OutputConnectionRepresentation::output_connection_t lineType = std::get<3>(externalConnectionLines_[i]);
        int lineColumn = std::get<2>(externalConnectionLines_[i]);
        verticalLinesInCurrentRow.push_back(std::make_tuple(lineColumn, false, false, lineType));
      }
    }

    // sort all found lines according to their column
    std::sort(verticalLinesInCurrentRow.begin(), verticalLinesInCurrentRow.end(), [](
      const std::tuple<int,bool,bool,solver_t::OutputConnectionRepresentation::output_connection_t> &a,
      const std::tuple<int,bool,bool,solver_t::OutputConnectionRepresentation::output_connection_t> &b)
    {
      return std::get<0>(a) < std::get<0>(b);
    });

    VLOG(1) << "line " << currentLineNo << ", currentLineNo " << currentLineNo << ", line=[" << line << "] has vertical lines: " << verticalLinesInCurrentRow;

    // fill current line with spaces such that it has the correct length
    const int requiredLineLength = VARIABLES_LINE_LENGTH + nColumnsForInternalConnectionLines + 5;  //< line length for structure and variables
    lineLength = StringUtility::stringLength(line);

    // determine fill character for horizontal line
    std::string fillCharacter = " ";
    if (lineStartsHere || lineEndsHere)
    {
      if (lineType == solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalReuse)
      {
        fillCharacter = "═";
      }
      else
      {
        fillCharacter = "─";
      }
    }
    else if (std::find(slotLineNosAll_.begin(), slotLineNosAll_.end(), currentLineNo) != slotLineNosAll_.end())
    {
      unconnectedSlotHere = true;
    }

    // if the line is yet too short, add fillCharacter until is has the correct length
    if (requiredLineLength > lineLength)
    {
      for (int i = 0; i < requiredLineLength - lineLength; i++)
      {
        // starting tip of arrow if the line starts or ends in the current row and has the according type
        if (i == 0 && lineStartsHere && (lineType == solver_t::OutputConnectionRepresentation::output_connection_t::ba
          || lineType == solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalCopy))
        {
          line += "<";
        }
        else if (i == 0 && lineEndsHere && (lineType == solver_t::OutputConnectionRepresentation::output_connection_t::ab
          || lineType == solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalCopy))
        {
          line += "<";
        }
        else
        {
          line += fillCharacter;
        }
      }
    }

    // add external connection lines
    // ---------
    // loop over the vertical lines in the current row
    int linesIndex = 0;
    if (!verticalLinesInCurrentRow.empty())
    {
      // loop over columns in the current row, current column is requiredLineLength + 2*i, i.e. it advances always by 2 characters
      for (int i = 0; i <= std::get<0>(verticalLinesInCurrentRow.back()); i++)
      {

        // determine events where line starts or ends
        bool lineCrosses = (lineStartsHere || lineEndsHere) && i < lineColumn;
        bool lineCurveStart = lineStartsHere && i == lineColumn;
        bool lineCurveEnd = lineEndsHere && i == lineColumn;

        // determine space character which is a horizontal line if the connection starts or ends here
        std::string space = " ";
        if (lineCrosses || ( (lineStartsHere || lineEndsHere) && i <= lineColumn))
        {
          if (lineType == solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalReuse)
          {
            space = "═";
          }
          else
          {
            space = "─";
          }
        }
        else if (unconnectedSlotHere && i == 0)
        {
          space = "x";
        }

        // if there is a vertical line passing at the current position
        if (i == std::get<0>(verticalLinesInCurrentRow[linesIndex]))
        {
          // depending on the type of the vertical line choose the right sign
          switch (std::get<3>(verticalLinesInCurrentRow[linesIndex]))
          {
          case solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalReuse:
            line += space;
            if (lineCrosses && !std::get<1>(verticalLinesInCurrentRow[linesIndex]) && !std::get<2>(verticalLinesInCurrentRow[linesIndex]))
            {
              line += "╬";
            }
            else if (lineCurveStart)
            {
              line += "╗";
            }
            else if (lineCurveEnd)
            {
              line += "╝";
            }
            else if (std::get<1>(verticalLinesInCurrentRow[linesIndex]))
            {
              line += "╦";
            }
            else if (std::get<2>(verticalLinesInCurrentRow[linesIndex]))
            {
              line += "╩";
            }
            else
            {
              line += "║";
            }
            break;
          case solver_t::OutputConnectionRepresentation::output_connection_t::ab:
          case solver_t::OutputConnectionRepresentation::output_connection_t::ba:
          case solver_t::OutputConnectionRepresentation::output_connection_t::bidirectionalCopy:
            line += space;
            if (lineCrosses && !std::get<1>(verticalLinesInCurrentRow[linesIndex]) && !std::get<2>(verticalLinesInCurrentRow[linesIndex]))
            {
              line += "┼";
            }
            else if (lineCurveStart)
            {
              line += "┐";
            }
            else if (lineCurveEnd)
            {
              line += "┘";
            }
            else if (std::get<1>(verticalLinesInCurrentRow[linesIndex]))
            {
              line += "┬";
            }
            else if (std::get<2>(verticalLinesInCurrentRow[linesIndex]))
            {
              line += "┴";
            }
            else
            {
              line += "│";
            }
            break;
          };
          linesIndex++;
        }
        else if (unconnectedSlotHere && i == 0)
        {
          line += std::string("x ");
        }
        else
        {
          // there is no vertical line passing through the current position, add two spaces
          line += space;
          line += space;
        }
      }
    }
    else if (unconnectedSlotHere)
    {
      line += std::string("x");
    }

    // add endline at the end of the current line
    line += "\n";

    // add the code line to the diagram
    diagram += line;
  }

  diagram += "connection types:\n";
  diagram += "  +··+   internal connection, no copy\n";
  diagram += "  ════   reuse variable, no copy\n";
  diagram += "  ───>   copy data in direction of arrow\n";

  //LOG(DEBUG) << "solver structure without connections:\n" << diagramWithoutConnectionLines.str() << "\n";
  LOG(DEBUG) << "solver structure with connections:\n" << diagram << "\n";

  return diagram;
}

//! produce the resulting file
void SolverStructureVisualizer::
writeDiagramFile(std::string filename)
{
  // only produce file on rank 0
  if (DihuContext::ownRankNoCommWorld() == 0)
  {
    std::string diagram = getDiagram();

    // output diagram to a text file
    std::ofstream file;
    OutputWriter::Generic::openFile(file, filename);

    if (!file.is_open())
    {
      LOG(FATAL) << "Could not write to file \"" << filename << "\".";
    }

    // output diagram text
    file << diagram;
    file.close();

    LOG(INFO) << "File \"" << filename << "\" written.";
  }

  LOG(DEBUG) << "SolverStructureVisualizer::clear data";

  // clear all values, this is required for unit tests where multiple different setups are created after each other
  solverRoot_ = nullptr;
  currentSolver_ = nullptr;
}
