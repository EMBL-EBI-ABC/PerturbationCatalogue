import {Component, ViewChild, AfterViewInit} from '@angular/core';
import {MatPaginator, MatPaginatorModule} from '@angular/material/paginator';
import {MatSort, MatSortModule, SortDirection} from '@angular/material/sort';
import {merge, of as observableOf} from 'rxjs';
import {catchError, map, startWith, switchMap} from 'rxjs/operators';
import {MatTableModule} from '@angular/material/table';
import {MatProgressSpinnerModule} from '@angular/material/progress-spinner';
import {MaveDBData, MaveDBFilters} from "../model/mavedb";
import {MatFormFieldModule} from "@angular/material/form-field";
import {MatInputModule} from "@angular/material/input";
import {MatCardModule} from "@angular/material/card";
import {MatListModule} from "@angular/material/list";
import {MatDividerModule} from "@angular/material/divider";
import {ApiService} from "../services/api.service";
import {RouterLink} from "@angular/router";

@Component({
  selector: 'app-data-portal',
  imports: [
    MatProgressSpinnerModule,
    MatTableModule,
    MatSortModule,
    MatPaginatorModule,
    MatFormFieldModule,
    MatInputModule,
    MatCardModule,
    MatListModule,
    MatDividerModule,
    RouterLink
  ],
  templateUrl: './data-portal.component.html',
  styleUrl: './data-portal.component.scss'
})
export class DataPortalComponent implements AfterViewInit {

  displayedColumns: string[] = [
    'urn',
    'sequenceType',
    'geneName',
    'geneCategory',
    'publicationYear',
    'numVariants'
  ];

  data: MaveDBData[] = [];
  aggregations: any;
  searchValue!: string;
  filters: MaveDBFilters = {
    sequenceType: [],
    geneCategory: [],
    publicationYear: []
  };

  resultsLength = 0;
  isLoadingResults = true;

  @ViewChild(MatPaginator) paginator!: MatPaginator;
  @ViewChild(MatSort) sort!: MatSort;

  constructor(private apiService: ApiService) {
  }

  ngAfterViewInit() {
    // If the user changes the sort order, reset back to the first page.
    this.sort.sortChange.subscribe(() => (this.paginator.pageIndex = 0));

    merge(this.sort.sortChange, this.paginator.page)
      .pipe(
        startWith({}),
        switchMap(() => {
          this.isLoadingResults = true;
          return this.apiService!.getAllMaveDBData(
            this.paginator.pageIndex,
            this.paginator.pageSize,
            this.sort.active,
            this.sort.direction,
            this.searchValue,
            this.filters,
          ).pipe(catchError(() => observableOf(null)));
        }),
        map(data => {
          // Flip flag to show that loading has finished.
          this.isLoadingResults = false;

          if (data === null) {
            return [];
          }

          // Only refresh the result length and aggregations if there is new data.
          this.resultsLength = data.total;
          this.aggregations = data.aggregations;
          return data.results;
        }),
      )
      .subscribe(data => (this.data = data));
  }

  search(event: Event) {
    this.searchValue = (event.target as HTMLInputElement).value.trim().toLowerCase();
    this.paginator.page.emit();
  }

  applyFilter(filterKey: string, filterValue: string) {
    // @ts-ignore
    const index = this.filters[filterKey].indexOf(filterValue);
    if (index > -1) {
      // @ts-ignore
      this.filters[filterKey].splice(index, 1)
    } else {
      // @ts-ignore
      this.filters[filterKey].push(filterValue);
    }
    this.paginator.page.emit();
  }

}
